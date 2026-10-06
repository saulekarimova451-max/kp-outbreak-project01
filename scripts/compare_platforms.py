"""
Compare Illumina (Snippy) and Nanopore (Clair3) variant calls on the same isolates.

Illumina is used as the reference standard ("truth") for SNPs because its per-base error rate is far
lower. For every isolate the script counts variants found by both platforms, by Illumina only and by
Nanopore only, separately for SNPs and indels, and records the homopolymer context of every
disagreement (Nanopore's typical error is a wrong length of a run like AAAAAA).

Usage (run by Snakemake rule 'compare_platforms'):
    python scripts/compare_platforms.py --ref resources/reference.fasta --isolates A,B,C \
        --illumina 'results/snippy/{iso}/snps.vcf' --ont 'results/ont/clair3/{iso}/merge_output.vcf.gz' \
        --summary report/tables/illumina_vs_nanopore.tsv \
        --discordant report/tables/illumina_vs_nanopore_discordant.tsv \
        --figure report/figures/illumina_vs_nanopore.png
"""
import argparse
import gzip
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

INDEL_WINDOW = 5  # bp: two callers may place the same indel at slightly different positions


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name:
        seqs[name] = "".join(chunks).upper()
    return seqs


def read_vcf(path, pass_only):
    """Return (snps, indels). SNPs as {(chrom, pos, alt)}; MNPs are split into single SNPs."""
    opener = gzip.open if str(path).endswith(".gz") else open
    snps, indels = set(), []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alts, _, flt = line.split("\t")[:7]
            if pass_only and flt not in ("PASS", "."):
                continue
            pos = int(pos)
            for alt in alts.split(","):
                if alt in ("*", "."):
                    continue
                if len(ref) == len(alt):
                    for i, (r, a) in enumerate(zip(ref, alt)):
                        if r != a:
                            snps.add((chrom, pos + i, a))
                else:
                    indels.append((chrom, pos, len(alt) - len(ref)))
    return snps, indels


def homopolymer(seq, pos):
    """Longest run of one base touching positions pos-1..pos+1 (1-based)."""
    best = 1
    for p in (pos - 1, pos, pos + 1):
        i = p - 1
        if i < 0 or i >= len(seq):
            continue
        base, left, right = seq[i], i, i
        while left > 0 and seq[left - 1] == base:
            left -= 1
        while right < len(seq) - 1 and seq[right + 1] == base:
            right += 1
        best = max(best, right - left + 1)
    return best


def match_indels(a, b):
    """Indels in a with a same-size indel in b within INDEL_WINDOW bp."""
    by_chrom = {}
    for c, p, d in b:
        by_chrom.setdefault(c, []).append((p, d))
    hit = []
    for c, p, d in a:
        hit.append(any(abs(p - q) <= INDEL_WINDOW and d == e for q, e in by_chrom.get(c, [])))
    return hit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--isolates", required=True)
    ap.add_argument("--illumina", required=True)
    ap.add_argument("--ont", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--discordant", required=True)
    ap.add_argument("--figure", required=True)
    args = ap.parse_args()

    ref = read_fasta(args.ref)
    summary, discordant = [], []
    for iso in args.isolates.split(","):
        ill_path, ont_path = Path(args.illumina.format(iso=iso)), Path(args.ont.format(iso=iso))
        if not (ill_path.exists() and ont_path.exists()):
            print(f"skip {iso}: missing input")
            continue
        ill_snps, ill_indels = read_vcf(ill_path, pass_only=False)  # Snippy writes only filtered calls
        ont_snps, ont_indels = read_vcf(ont_path, pass_only=True)

        shared = ill_snps & ont_snps
        ill_only, ont_only = ill_snps - ont_snps, ont_snps - ill_snps
        ill_hit = match_indels(ill_indels, ont_indels)
        ont_hit = match_indels(ont_indels, ill_indels)

        for platform, items in (("Illumina only", ill_only), ("Nanopore only", ont_only)):
            for c, p, a in sorted(items):
                discordant.append({"isolate": iso, "type": "SNP", "found_by": platform, "contig": c,
                                   "pos": p, "change": f"{ref[c][p - 1]}>{a}",
                                   "homopolymer_len": homopolymer(ref[c], p)})
        for platform, items, hits in (("Illumina only", ill_indels, ill_hit),
                                      ("Nanopore only", ont_indels, ont_hit)):
            for (c, p, d), h in zip(items, hits):
                if not h:
                    discordant.append({"isolate": iso, "type": "indel", "found_by": platform, "contig": c,
                                       "pos": p, "change": f"{d:+d} bp",
                                       "homopolymer_len": homopolymer(ref[c], p + 1)})

        n_shared_ind = sum(ont_hit)
        summary.append({
            "isolate": iso,
            "snps_illumina": len(ill_snps), "snps_nanopore": len(ont_snps), "snps_shared": len(shared),
            "snps_illumina_only": len(ill_only), "snps_nanopore_only": len(ont_only),
            "snp_precision": round(len(shared) / len(ont_snps), 3) if ont_snps else None,
            "snp_recall": round(len(shared) / len(ill_snps), 3) if ill_snps else None,
            "indels_illumina": len(ill_indels), "indels_nanopore": len(ont_indels),
            "indels_shared": n_shared_ind,
            "indels_nanopore_only": len(ont_indels) - n_shared_ind,
        })

    if not summary:
        raise SystemExit("No isolate has both an Illumina and a Nanopore VCF yet.")
    s = pd.DataFrame(summary)
    d = pd.DataFrame(discordant, columns=["isolate", "type", "found_by", "contig", "pos", "change",
                                          "homopolymer_len"])
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    s.to_csv(args.summary, sep="\t", index=False)
    d.to_csv(args.discordant, sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), gridspec_kw={"width_ratios": [3, 3, 2]})
    x = range(len(s))
    for ax, kind in zip(axes[:2], ("snps", "indels")):
        shared = s[f"{kind}_shared"]
        ill_only = s[f"{kind}_illumina"] - shared
        ont_only = s[f"{kind}_nanopore"] - shared
        ax.bar(x, shared, color="#2e7d32", label="both")
        ax.bar(x, ill_only, bottom=shared, color="#1565c0", label="Illumina only")
        ax.bar(x, ont_only, bottom=shared + ill_only, color="#ef6c00", label="Nanopore only")
        ax.set_xticks(list(x), s["isolate"], rotation=45, ha="right", fontsize=8)
        ax.set_title("SNPs" if kind == "snps" else "Indels", fontsize=11)
        ax.set_ylabel("variants vs reference")
    axes[0].legend(fontsize=8, frameon=False)
    if not d.empty:
        hp = d.groupby(["type", "homopolymer_len"]).size().unstack(0, fill_value=0)
        hp.plot(kind="bar", ax=axes[2], color={"SNP": "#6a1b9a", "indel": "#ef6c00"}, width=0.8)
        axes[2].set_xlabel("homopolymer length at the disagreement")
        axes[2].set_ylabel("disagreements")
        axes[2].set_title("Where the platforms disagree", fontsize=11)
        axes[2].legend(fontsize=8, frameon=False)
    fig.suptitle("Illumina (Snippy) vs Nanopore (Clair3) variant calls, same isolates, same reference", fontsize=11)
    fig.tight_layout()
    fig.savefig(args.figure, dpi=200)

    print(s.to_string(index=False))
    if not d.empty:
        print("\nDisagreements by type and homopolymer length:")
        print(d.groupby(["type", "found_by", "homopolymer_len"]).size().to_string())


if __name__ == "__main__":
    main()
