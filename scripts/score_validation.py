"""
Score Snippy's calls on simulated reads against the known truth: TP, FP, FN, precision, recall, F1.

Usage (run by Snakemake rule 'score_validation'):
    python scripts/score_validation.py results/validation/truth.tsv results/validation/snippy/snps.vcf \
        benchmarks/validation_snippy.tsv report/tables/validation_snippy.tsv
"""
import sys

import pandas as pd

from compare_platforms import read_vcf

truth_path, vcf_path, bench_path, out_path = sys.argv[1:5]
truth = pd.read_csv(truth_path, sep="\t")
true_snps = {(r.contig, int(r.pos), r.alt) for r in truth.itertuples()}
called, _ = read_vcf(vcf_path, pass_only=False)
tp, fp, fn = len(true_snps & called), len(called - true_snps), len(true_snps - called)
precision = tp / (tp + fp) if tp + fp else 0.0
recall = tp / (tp + fn) if tp + fn else 0.0
f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
runtime = pd.read_csv(bench_path, sep="\t")["s"].iloc[0]
out = pd.DataFrame([{"simulated_snps": len(true_snps), "called_snps": len(called), "true_positives": tp,
                     "false_positives": fp, "false_negatives": fn, "precision": round(precision, 4),
                     "recall": round(recall, 4), "f1": round(f1, 4), "snippy_runtime_s": round(runtime, 1)}])
out.to_csv(out_path, sep="\t", index=False)
print(out.to_string(index=False))
missed = sorted(true_snps - called)
if missed:
    print("\nMissed SNPs (contig, pos, alt):", missed[:20])
