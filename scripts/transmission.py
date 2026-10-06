"""
Apply the pre-registered SNP thresholds to the SNP distance matrix and draw the report figures.

Usage (run by Snakemake rule 'transmission', or by hand from the project folder):
    python scripts/transmission.py results/core/snp_matrix_clean.tsv results/tree/outbreak.treefile \
        config/samples.tsv config/config.yaml

Writes (committed to the repository, they are small):
    report/tables/snp_pairs.tsv        every pair: SNPs, days between samples, call (linked / grey zone / unrelated)
    report/tables/clusters.tsv         isolates grouped by single linkage at <= snp_linked SNPs
    report/figures/snp_heatmap.png     SNP distance matrix, isolates in date order
    report/figures/timeline.png        isolates on a time axis, each joined to its genetically closest earlier isolate
    report/figures/tree.png            maximum-likelihood tree (midpoint rooted) with dates and sample type
"""
import sys
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
import yaml
from Bio import Phylo


def load(matrix_path, samples_path, config_path):
    m = pd.read_csv(matrix_path, sep="\t", index_col=0)
    m.columns = m.columns.str.strip()
    m.index = m.index.str.strip()
    m = m.drop(index="Reference", columns="Reference", errors="ignore")
    s = pd.read_csv(samples_path, sep="\t", dtype=str).set_index("isolate")
    s["date"] = pd.to_datetime(s["collection_date"], errors="coerce")
    cfg = yaml.safe_load(open(config_path))["transmission"]
    order = s.loc[m.index].sort_values("date").index.tolist()
    return m.loc[order, order], s, cfg["snp_linked"], cfg["snp_unrelated"]


def call(snps, linked, unrelated):
    if snps <= linked:
        return "linked"
    return "grey zone" if snps <= unrelated else "unrelated"


def clusters(m, linked):
    """Single-linkage clustering: isolates joined if any chain of pairs is <= linked SNPs."""
    parent = {i: i for i in m.index}

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for a, b in combinations(m.index, 2):
        if m.loc[a, b] <= linked:
            parent[find(a)] = find(b)
    roots = {r: n for n, r in enumerate(sorted({find(i) for i in m.index}, key=lambda r: list(m.index).index(r)), 1)}
    return {i: roots[find(i)] for i in m.index}


def heatmap(m, s, linked, out):
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    im = ax.imshow(m.values, cmap="Blues", vmin=0, vmax=max(linked, int(m.values.max())))
    labels = [f"{i}\n{s.loc[i, 'date']:%d %b}" for i in m.index]
    ax.set_xticks(range(len(m)), labels, rotation=90, fontsize=8)
    ax.set_yticks(range(len(m)), labels, fontsize=8)
    for r in range(len(m)):
        for c in range(len(m)):
            v = int(m.iat[r, c])
            ax.text(c, r, v, ha="center", va="center", fontsize=9,
                    color="white" if v > 0.6 * max(linked, m.values.max()) else "black")
    fig.colorbar(im, ax=ax, shrink=0.8, label="SNP differences")
    ax.set_title(f"Pairwise core-genome SNP distances (linked if <= {linked})", fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def timeline(m, s, out):
    """Two panels with a broken time axis: the dense April-May cluster and the late isolate(s)."""
    order = list(m.index)  # already sorted by date
    y = {iso: n for n, iso in enumerate(order)}
    dates = s.loc[order, "date"]
    gap = dates.diff().dt.days.fillna(0)
    split = gap.idxmax() if gap.max() > 21 else None  # first isolate after a gap of > 3 weeks
    early = order if split is None else order[: order.index(split)]
    late = [] if split is None else order[order.index(split):]

    pad = pd.Timedelta(days=2)
    label_room = pd.Timedelta(days=3)  # space for the isolate name right of the last point
    panels = [(early, dates[early].min() - pad, dates[early].max() + label_room)]
    if late:
        panels.append((late, dates[late].min() - pad, dates[late].max() + label_room))
    widths = [max((hi - lo).days, 4) for _, lo, hi in panels]
    fig, axes = plt.subplots(1, len(panels), figsize=(10, 5), sharey=True,
                             gridspec_kw={"width_ratios": widths, "wspace": 0.05})
    axes = list(axes) if len(panels) > 1 else [axes]
    colors = {"urine": "#c0392b", "rectal swab": "#2471a3"}

    for ax, (members, lo, hi) in zip(axes, panels):
        ax.set_xlim(lo, hi)
        for iso in members:
            src = str(s.loc[iso, "isolation_source"]).lower()
            ax.scatter(dates[iso], y[iso], s=60, color=colors.get(src, "0.3"), zorder=3)
            ax.text(dates[iso] + pd.Timedelta(hours=10), y[iso], f" {iso}", fontsize=8, va="center")
        if members is late:
            ax.set_xticks(sorted(set(dates[members])))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
        ax.tick_params(axis="x", labelrotation=45, labelsize=8)
        ax.set_yticks([])

    for n, iso in enumerate(order[1:], 1):
        nearest = min(order[:n], key=lambda e: m.loc[iso, e])
        d = int(m.loc[iso, nearest])
        same_panel = (iso in early) == (nearest in early)
        if same_panel:
            ax = axes[0] if iso in early else axes[-1]
            ax.plot([dates[nearest], dates[iso]], [y[nearest], y[iso]], color="0.6", lw=1, zorder=1)
            ax.text(dates[iso] - pd.Timedelta(hours=14), y[iso] - 0.35, f"{d} SNP" + ("s" if d != 1 else ""),
                    fontsize=7, color="0.25", ha="right")
        else:
            axes[-1].text(dates[iso], y[iso] - 0.45, f"{d} SNPs from {nearest}", fontsize=7,
                          color="0.25", ha="center")

    if len(axes) > 1:
        axes[0].spines["right"].set_visible(False)
        axes[1].spines["left"].set_visible(False)
        for ax, x in ((axes[0], 1), (axes[1], 0)):
            ax.plot([x, x], [0, 1], transform=ax.transAxes, color="white", lw=3, clip_on=False)
            ax.text(x, 0, "//", transform=ax.transAxes, ha="center", va="center", fontsize=10)
    for src, col in colors.items():
        axes[0].scatter([], [], color=col, label=src)
    axes[0].legend(loc="upper left", fontsize=8, frameon=False)
    axes[0].set_ylim(-0.8, len(order) - 0.4)
    fig.suptitle("Isolates by collection date; grey line = SNPs to the genetically closest earlier isolate",
                 fontsize=10)
    fig.subplots_adjust(bottom=0.18, top=0.9, left=0.04, right=0.97)
    fig.savefig(out, dpi=200)
    plt.close(fig)


def tree_plot(tree_path, s, out):
    tree = Phylo.read(tree_path, "newick")
    # 'Reference' is the assembly of isolate 0326576 itself; it adds no information to the tree
    for tip in tree.get_terminals():
        if tip.name == "Reference":
            tree.prune(tip)
    tree.root_at_midpoint()
    tree.ladderize()
    for tip in tree.get_terminals():
        if tip.name in s.index:
            tip.name = f"{tip.name}  ({s.loc[tip.name, 'date']:%d %b}, {str(s.loc[tip.name, 'isolation_source']).lower()})"
    fig, ax = plt.subplots(figsize=(9, 5))
    Phylo.draw(tree, axes=ax, do_show=False,
               label_func=lambda c: c.name if c.is_terminal() else None,
               branch_labels=lambda c: (f"{c.confidence:.0f}" if not c.is_terminal() and c.confidence is not None
                                        and c.confidence >= 70 else None))
    depth = max(tree.depths().values())
    ax.set_xlim(-0.02 * depth, depth * 1.6)  # room for tip labels
    ax.set_xlabel("Substitutions per variable site")
    ax.set_ylabel("")
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.set_title("Maximum-likelihood tree (IQ-TREE, midpoint rooted); numbers = bootstrap support >= 70",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main():
    matrix_path, tree_path, samples_path, config_path = sys.argv[1:5]
    tables, figures = Path("report/tables"), Path("report/figures")
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    m, s, linked, unrelated = load(matrix_path, samples_path, config_path)

    rows = []
    for a, b in combinations(m.index, 2):
        snps = int(m.loc[a, b])
        days = abs((s.loc[b, "date"] - s.loc[a, "date"]).days)
        rows.append({"isolate_a": a, "isolate_b": b, "snps": snps, "days_apart": days,
                     "call": call(snps, linked, unrelated)})
    pairs = pd.DataFrame(rows).sort_values(["snps", "days_apart"])
    pairs.to_csv(tables / "snp_pairs.tsv", sep="\t", index=False)

    cl = clusters(m, linked)
    pd.DataFrame({"isolate": list(cl), "cluster": list(cl.values()),
                  "collection_date": [s.loc[i, "collection_date"] for i in cl],
                  "isolation_source": [s.loc[i, "isolation_source"] for i in cl]}
                 ).to_csv(tables / "clusters.tsv", sep="\t", index=False)

    heatmap(m, s, linked, figures / "snp_heatmap.png")
    timeline(m, s, figures / "timeline.png")
    tree_plot(tree_path, s, figures / "tree.png")

    print(f"{len(m)} isolates, {len(pairs)} pairs, {len(set(cl.values()))} cluster(s) at <= {linked} SNPs")
    print(pairs["call"].value_counts().to_string())
    print(f"SNP distance: min {pairs.snps.min()}, median {pairs.snps.median():.0f}, max {pairs.snps.max()}")


if __name__ == "__main__":
    main()
