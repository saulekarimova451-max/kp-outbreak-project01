"""
Make a ground-truth test for SNP calling: copy the reference and change N random bases.

The simulated isolate then differs from the reference at exactly these N known positions, so the
variant caller's precision and recall can be measured instead of assumed (rubric criterion 6).

Usage (run by Snakemake rule 'mutate_reference'):
    python scripts/mutate_reference.py resources/reference.fasta results/validation/mutated.fasta \
        results/validation/truth.tsv 100 42
"""
import random
import sys

from compare_platforms import read_fasta

ref_path, out_fa, out_truth, n, seed = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
rng = random.Random(seed)
seqs = read_fasta(ref_path)
# only contigs >= 20 kb, and not within 500 bp of a contig end (reads cannot cover ends evenly)
eligible = [(c, len(s)) for c, s in seqs.items() if len(s) >= 20000]
weights = [L for _, L in eligible]
mutable = {c: list(s) for c, s in seqs.items()}
truth, used = [], set()
while len(truth) < n:
    c, L = rng.choices(eligible, weights=weights)[0]
    p = rng.randint(501, L - 500)
    if (c, p) in used or any(abs(p - q) < 50 for cc, q in used if cc == c):
        continue  # keep SNPs >= 50 bp apart so each is called as a single SNP
    ref = mutable[c][p - 1]
    if ref not in "ACGT":
        continue
    alt = rng.choice([b for b in "ACGT" if b != ref])
    mutable[c][p - 1] = alt
    used.add((c, p))
    truth.append((c, p, ref, alt))

with open(out_fa, "w") as f:
    for c, s in mutable.items():
        f.write(f">{c}\n")
        seq = "".join(s)
        for i in range(0, len(seq), 80):
            f.write(seq[i:i + 80] + "\n")
with open(out_truth, "w") as f:
    f.write("contig\tpos\tref\talt\n")
    for c, p, r, a in sorted(truth):
        f.write(f"{c}\t{p}\t{r}\t{a}\n")
print(f"{n} SNPs written to {out_truth}")
