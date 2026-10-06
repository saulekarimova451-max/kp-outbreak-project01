"""
Summarise Snakemake benchmark files: runtime per step, and a scale estimate for 100 and 1,000 isolates.

Usage (run by Snakemake rule 'runtime_summary'):
    python scripts/runtime_summary.py benchmarks report/tables/runtime.tsv report/figures/runtime.png 8 4
    (arguments: benchmark folder, output table, output figure, number of isolates, cores used)
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

bench_dir, out_table, out_fig = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
n_isolates, cores = int(sys.argv[4]), int(sys.argv[5])

rows = []
for f in sorted(bench_dir.rglob("*.tsv")):
    rule = f.parent.name if f.parent != bench_dir else f.stem
    try:
        b = pd.read_csv(f, sep="\t")
    except Exception:
        continue
    if b.empty or "s" not in b:
        continue
    rss = pd.to_numeric(b.get("max_rss", pd.Series(["NA"])), errors="coerce").iloc[0]
    rows.append({"rule": rule, "job": f.stem, "seconds": float(b["s"].iloc[0]), "max_rss_mb": rss})
d = pd.DataFrame(rows)
per_rule = (d.groupby("rule")
            .agg(jobs=("job", "count"), mean_s=("seconds", "mean"), total_min=("seconds", lambda s: s.sum() / 60),
                 max_rss_mb=("max_rss_mb", "max"))
            .round(1).sort_values("total_min", ascending=False))
# Per-isolate steps scale linearly with the number of isolates; one-off steps (gubbins, iqtree, databases) do not.
per_isolate = per_rule[per_rule.jobs >= n_isolates]
minutes_per_isolate = per_isolate["total_min"].sum() / n_isolates
per_rule.to_csv(out_table, sep="\t")

fig, ax = plt.subplots(figsize=(7, 0.35 * len(per_rule) + 1.2))
ax.barh(per_rule.index[::-1], per_rule["total_min"][::-1], color="#1565c0")
ax.set_xlabel("total wall-clock minutes (all isolates)")
ax.set_title(f"Runtime per workflow step ({n_isolates} isolates, {cores} cores)", fontsize=10)
fig.tight_layout()
fig.savefig(out_fig, dpi=200)

print(per_rule.to_string())
print(f"\nPer-isolate steps: {minutes_per_isolate:.1f} min per isolate on {cores} cores")
for n in (100, 1000):
    hours = minutes_per_isolate * n / 60
    print(f"Estimate for {n} isolates: {hours:.1f} h on this laptop; "
          f"{hours / 25:.1f} h on a 100-core cluster (25 jobs of 4 cores in parallel)")
