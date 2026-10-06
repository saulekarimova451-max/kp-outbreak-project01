"""
Would a Nanopore-only analysis reach the same outbreak conclusion as Illumina?

1. Pairwise SNP distances between isolates, computed the same way from Illumina (Snippy) and from
   Nanopore (Clair3) variant calls: distance = SNPs present in one isolate but not the other.
2. Profile of the Nanopore-only SNPs: how many sit in bacterial methylation motifs (GATC = Dam,
   CCWGG = Dcm), compared with the share of the genome covered by those motifs, and how many recur
   in every isolate (a systematic error, which cancels out when isolates are compared).

Usage (run by Snakemake rule 'nanopore_distances'):
    python scripts/nanopore_distances.py --ref resources/reference.fasta --isolates A,B,... \
        --illumina 'results/snippy/{iso}/snps.vcf' --ont 'results/ont/clair3/{iso}/merge_output.vcf.gz' \
        --matrix report/tables/snp_matrix_nanopore.tsv --profile report/tables/nanopore_only_snp_profile.tsv \
        --figure report/figures/illumina_vs_nanopore_distances.png
"""
import argparse
import re
from collections import Counter
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from compare_platforms import read_fasta, read_vcf

MOTIFS = {"GATC (Dam)": "GATC", "CCWGG (Dcm)": "CC[AT]GG"}


def motif_mask(seqs):
    """Positions (contig, 1-based) covered by each methylation motif."""
    masks = {}
    for name, pattern in MOTIFS.items():
        covered = set()
        rx = re.compile(f"(?=({pattern}))")
        for c, s in seqs.items():
            for m in rx.finditer(s):
                start = m.start() + 1
                covered.update((c, p) for p in range(start, start + len(m.group(1))))
        masks[name] = covered
    return masks


def distances(calls, isolates):
    m = pd.DataFrame(0, index=isolates, columns=isolates)
    for a, b in combinations(isolates, 2):
        d = len(calls[a] ^ calls[b])
        m.loc[a, b] = m.loc[b, a] = d
    return m


def main():
    ap = argparse.ArgumentParser()
    for arg in ("ref", "isolates", "illumina", "ont", "matrix", "profile", "figure"):
        ap.add_argument(f"--{arg}", required=True)
    args = ap.parse_args()

    ref = read_fasta(args.ref)
    genome_len = sum(len(s) for s in ref.values())
    isolates = [i for i in args.isolates.split(",")
                if Path(args.illumina.format(iso=i)).exists() and Path(args.ont.format(iso=i)).exists()]
    ill = {i: read_vcf(args.illumina.format(iso=i), pass_only=False)[0] for i in isolates}
    ont = {i: read_vcf(args.ont.format(iso=i), pass_only=True)[0] for i in isolates}

    m_ill, m_ont = distances(ill, isolates), distances(ont, isolates)
    out = m_ont.copy()
    out.to_csv(args.matrix, sep="\t")

    # Nanopore-only SNPs and their context
    masks = motif_mask(ref)
    ont_only = {i: ont[i] - ill[i] for i in isolates}
    recurrence = Counter(snp for i in isolates for snp in ont_only[i])
    rows = []
    for i in isolates:
        snps = ont_only[i]
        row = {"isolate": i, "nanopore_only_snps": len(snps),
               "in_all_isolates": sum(recurrence[s] == len(isolates) for s in snps),
               "in_one_isolate_only": sum(recurrence[s] == 1 for s in snps)}
        for name, covered in masks.items():
            row[f"in_{name.split()[0]}"] = sum((c, p) in covered for c, p, _ in snps)
        rows.append(row)
    prof = pd.DataFrame(rows)
    background = {f"in_{n.split()[0]}": len(cov) / genome_len for n, cov in masks.items()}
    total = prof["nanopore_only_snps"].sum()
    summary = {"isolate": "ALL", "nanopore_only_snps": total,
               "in_all_isolates": prof["in_all_isolates"].sum(),
               "in_one_isolate_only": prof["in_one_isolate_only"].sum()}
    for k in background:
        summary[k] = prof[k].sum()
    prof = pd.concat([prof, pd.DataFrame([summary])], ignore_index=True)
    for k, bg in background.items():
        prof[f"{k}_pct"] = (100 * prof[k] / prof["nanopore_only_snps"]).round(1)
        prof[f"{k}_genome_background_pct"] = round(100 * bg, 1)
    prof.to_csv(args.profile, sep="\t", index=False)

    # Figure: Illumina vs Nanopore pairwise distances
    pairs = [(a, b) for a, b in combinations(isolates, 2)]
    x = [m_ill.loc[a, b] for a, b in pairs]
    y = [m_ont.loc[a, b] for a, b in pairs]
    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.scatter(x, y, color="#1565c0", alpha=0.8)
    lim = max(max(x), max(y)) * 1.1 + 1
    ax.plot([0, lim], [0, lim], color="0.6", ls="--", lw=1, label="same distance")
    ax.axhline(21, color="#c62828", lw=1, label="21-SNP threshold")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("SNP distance from Illumina calls")
    ax.set_ylabel("SNP distance from Nanopore calls")
    ax.set_title(f"Every isolate pair ({len(pairs)}), same isolates, both platforms", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(args.figure, dpi=200)

    print("Nanopore pairwise SNP distances:\n" + m_ont.to_string())
    print(f"\nPairs: Illumina max {max(x)}, Nanopore max {max(y)}; "
          f"pairs above 21 SNPs with Nanopore: {sum(v > 21 for v in y)} of {len(y)}")
    print("\nNanopore-only SNP profile:\n" + prof.to_string(index=False))


if __name__ == "__main__":
    main()
