"""
Turn BLAST hits against CARD into a resistance-gene table, a predicted antibiogram, and a comparison
with NCBI Pathogen Detection's AMR genotypes.

Filtering (justified in the report):
- identity >= 90%: below this a hit is usually a different gene variant or family;
- coverage of the reference gene >= 80%: below this the gene is likely broken or split across contigs;
- overlapping hits on the same contig region: only the best (highest bitscore) is kept.

Drug-class prediction uses explicit gene-name rules (RULES below), not CARD's broad family-level
drug classes, because those would e.g. call every SHV beta-lactamase a carbapenemase. Efflux pumps and
regulators are listed but not used for the prediction (they are intrinsic to K. pneumoniae).

Usage (run by Snakemake rule 'amr_summary'):
    python scripts/amr_summary.py --blast results/amr/blast --aro resources/card/aro_index.tsv \
        --samples config/samples.tsv --ncbi report/tables/pathogen_detection_isolates.tsv \
        --flye results/ont/flye --outdir report
"""
import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

MIN_ID, MIN_COV = 90.0, 80.0

# (drug class, regex on the CARD model name). Order = column order in the antibiogram.
RULES = [
    ("Carbapenems", r"^(KPC|NDM|OXA-(48|181|232|244|162|204)|VIM|IMP|GES-(2|4|5|14))\b"),
    ("3rd-gen cephalosporins", r"^(KPC|NDM|VIM|IMP|OXA-(48|181|232)|CTX-M|SHV-(2|5|12|27|31)\b|CMY|DHA|ACC|FOX|GES|PER|VEB)"),
    ("Aminoglycosides", r"^(AAC\(|ANT\(|APH\(|aad|rmt|armA|npmA)"),
    ("Fluoroquinolones", r"^(qnr|QnrB|QnrS|QnrA|AAC\(6'\)-Ib-cr)"),
    ("Sulfonamides", r"^sul[123]"),
    ("Trimethoprim", r"^dfr"),
    ("Tetracyclines", r"^tet\("),
    ("Phenicols", r"^(cat[AB]|floR|cmlA)"),
    ("Macrolides", r"^(mph|erm|msr|ere)"),
    ("Rifampicin", r"^arr"),
    ("Colistin", r"^MCR-"),
]
INTRINSIC = r"(efflux|regulator)"
# Chromosomal (intrinsic) K. pneumoniae genes: present in every isolate, not acquired, not used for prediction.
INTRINSIC_NAMES = (r"^(arnT|eptB|ompA|OmpK|Klebsiella pneumoniae|Kpn[EFGH]|mdt|fosA|oqx|emrD|acr|marA|ramA|soxS|"
                   r"tolC|msbA|H-NS|baeR|cpxA|golS|LEN-|OKP-|SHV-(1|11|26|28|182)\b)")
# Same gene, different name in CARD vs NCBI AMRFinderPlus
SYNONYMS = {"brpmbl": "ble", "mrxa": "mrx"}
NCBI_INTRINSIC = {"fosa", "oqxa", "oqxb", "emrd", "shv"}


def card_name(sseqid, aro_index):
    """CARD headers look like gb|AF0000.1|+|0-861|ARO:3002312|KPC-2 [Klebsiella pneumoniae].
    The database build replaces '|' with '__' and spaces with '_' so BLAST keeps the whole header.
    The gene name is taken from CARD's own index (ARO Name), not from the header."""
    parts = re.split(r"__|\|", sseqid)
    aro = next((p for p in parts if p.startswith("ARO:")), "")
    if aro in aro_index.index:
        return aro, aro_index.loc[aro, "ARO Name"]
    return aro, parts[-1].split("_[")[0]


def best_hits(df):
    """Greedy: keep the highest-bitscore hit, drop others overlapping it by > 50% on the same contig."""
    keep = []
    for _, h in df.sort_values("bitscore", ascending=False).iterrows():
        lo, hi = sorted((h.qstart, h.qend))
        clash = False
        for k in keep:
            if k["qseqid"] != h.qseqid:
                continue
            klo, khi = sorted((k["qstart"], k["qend"]))
            overlap = min(hi, khi) - max(lo, klo)
            if overlap > 0.5 * min(hi - lo, khi - klo):
                clash = True
                break
        if not clash:
            keep.append(h)
    return pd.DataFrame(keep)


def norm(gene):
    g = gene.lower().replace("bla", "")
    g = re.sub(r"[^a-z0-9]", "", g)
    return SYNONYMS.get(g, g)


def gene_family(gene):
    """Gene family = name without the allele number: blaNDM-5 -> ndm, aac(6')-Ib9 -> aac6ib, sul1 -> sul."""
    g = norm(gene)
    g = re.sub(r"(?<=[a-z])\d+$", "", g) if not g.startswith(("aac", "aph", "ant")) else re.sub(r"\d+$", "", g)
    return SYNONYMS.get(g, g)


def main():
    ap = argparse.ArgumentParser()
    for arg in ("blast", "aro", "samples", "ncbi", "flye", "outdir"):
        ap.add_argument(f"--{arg}", required=True)
    args = ap.parse_args()
    tables, figures = Path(args.outdir, "tables"), Path(args.outdir, "figures")
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    aro = pd.read_csv(args.aro, sep="\t", dtype=str)
    aro = aro.set_index("ARO Accession")
    cols = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "slen", "evalue", "bitscore"]

    genes = []
    for f in sorted(Path(args.blast).glob("*.tsv")):
        sample = f.stem
        df = pd.read_csv(f, sep="\t", names=cols)
        if df.empty:
            continue
        df["coverage"] = 100 * df["length"] / df["slen"]
        n_raw = len(df)
        df = df[(df.pident >= MIN_ID) & (df.coverage >= MIN_COV)]
        n_pass = len(df)
        df = best_hits(df) if not df.empty else df
        print(f"{sample}: {n_raw} raw hits, {n_pass} pass identity/coverage, {len(df)} after overlap filter")
        info = None
        if sample.startswith("flye-"):
            p = Path(args.flye, sample[len("flye-"):], "assembly_info.txt")
            if p.exists():
                info = pd.read_csv(p, sep="\t").set_index("#seq_name")
        for _, h in df.iterrows():
            a, name = card_name(h.sseqid, aro)
            mech = aro.loc[a, "Resistance Mechanism"] if a in aro.index else ""
            family = aro.loc[a, "AMR Gene Family"] if a in aro.index else ""
            row = {"sample": sample, "contig": h.qseqid, "gene": name, "aro": a,
                   "identity": round(h.pident, 1), "coverage": round(h.coverage, 1),
                   "gene_family": family, "mechanism": mech,
                   "intrinsic_or_efflux": bool(re.search(INTRINSIC, f"{mech} {family}", re.I)
                                               or re.search(INTRINSIC_NAMES, name, re.I))}
            if info is not None and h.qseqid in info.index:
                row["contig_length"] = int(info.loc[h.qseqid, "length"])
                row["contig_circular"] = info.loc[h.qseqid, "circ."]
            genes.append(row)

    g = pd.DataFrame(genes)
    g.to_csv(tables / "amr_genes.tsv", sep="\t", index=False)

    # Predicted antibiogram (acquired genes only)
    acquired = g[~g.intrinsic_or_efflux]
    samples = list(dict.fromkeys(g["sample"]))
    abg = []
    for s in samples:
        row = {"sample": s}
        names = acquired.loc[acquired["sample"] == s, "gene"].tolist()
        for drug, rx in RULES:
            hits = sorted({n for n in names if re.search(rx, n)})
            row[drug] = "R (" + ", ".join(hits) + ")" if hits else "S*"
        abg.append(row)
    abg = pd.DataFrame(abg)
    abg.to_csv(tables / "antibiogram_predicted.tsv", sep="\t", index=False)

    # Comparison with NCBI AMRFinderPlus genotypes (Illumina assemblies only)
    smp = pd.read_csv(args.samples, sep="\t", dtype=str).set_index("illumina_biosample")
    ncbi = pd.read_csv(args.ncbi, sep="\t", dtype=str)
    comp = []
    for _, r in ncbi.iterrows():
        iso = smp.loc[r["BioSample"], "isolate"] if r["BioSample"] in smp.index else None
        if iso is None or iso not in samples:
            continue
        ncbi_names = [x.split("=")[0] for x in str(r["AMR genotypes"]).split(",") if x and not x.endswith("=POINT")]
        theirs = {norm(x) for x in ncbi_names if gene_family(x) not in NCBI_INTRINSIC}
        ours_names = acquired.loc[acquired["sample"] == iso, "gene"].tolist()
        ours = {norm(n) for n in ours_names}
        f_theirs = {gene_family(x) for x in ncbi_names if gene_family(x) not in NCBI_INTRINSIC}
        f_ours = {gene_family(n) for n in ours_names}
        comp.append({"isolate": iso,
                     "families_ncbi": len(f_theirs), "families_ours": len(f_ours),
                     "families_both": len(f_theirs & f_ours),
                     "families_ncbi_only": ", ".join(sorted(f_theirs - f_ours)),
                     "families_ours_only": ", ".join(sorted(f_ours - f_theirs)),
                     "alleles_both": len(theirs & ours),
                     "alleles_ncbi_only": ", ".join(sorted(theirs - ours)),
                     "alleles_ours_only": ", ".join(sorted(ours - theirs)),
                     "ncbi_point_mutations": ", ".join(x.split("=")[0] for x in str(r["AMR genotypes"]).split(",")
                                                      if x.endswith("=POINT"))})
    comp = pd.DataFrame(comp)
    comp.to_csv(tables / "amr_vs_ncbi.tsv", sep="\t", index=False)

    # Figure: presence/absence of acquired genes
    if not acquired.empty:
        mat = (acquired.assign(v=1).pivot_table(index="gene", columns="sample", values="v", aggfunc="max")
               .fillna(0).reindex(columns=samples, fill_value=0))
        carb = [n for n in mat.index if re.search(RULES[0][1], n)]
        mat = mat.loc[carb + [n for n in mat.index if n not in carb]]
        fig, ax = plt.subplots(figsize=(1.0 + 0.55 * len(samples), 0.6 + 0.28 * len(mat)))
        ax.imshow(mat.values, cmap="Greens", vmin=0, vmax=1.4, aspect="auto")
        ax.set_xticks(range(len(samples)), samples, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(len(mat)), mat.index, fontsize=8)
        for t in ax.get_yticklabels():
            if t.get_text() in carb:
                t.set_color("#c62828")
                t.set_fontweight("bold")
        ax.set_title("Acquired resistance genes (BLAST vs CARD; carbapenemases in red)", fontsize=10)
        fig.tight_layout()
        fig.savefig(figures / "amr_genes.png", dpi=200)

    print("\nPredicted antibiogram:\n" + abg.to_string(index=False))
    if len(comp):
        print("\nGene families vs NCBI AMRFinderPlus (acquired genes):")
        print(comp[["isolate", "families_ncbi", "families_ours", "families_both", "families_ncbi_only",
                    "families_ours_only"]].to_string(index=False))
    carb_rows = g[g.gene.apply(lambda n: bool(re.search(RULES[0][1], n)))]
    if not carb_rows.empty:
        print("\nCarbapenemase genes and their contigs:")
        print(carb_rows.drop(columns=["aro", "gene_family", "mechanism", "intrinsic_or_efflux"]).to_string(index=False))


if __name__ == "__main__":
    main()
