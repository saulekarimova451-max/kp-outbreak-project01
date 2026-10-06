# Phase 5.5-5.6: apply the pre-registered thresholds and draw the report figures.


rule transmission:
    """Pairs table, clusters at <= snp_linked SNPs, heatmap, timeline and tree figures."""
    input:
        matrix="results/core/snp_matrix_clean.tsv",
        tree="results/tree/outbreak.treefile",
        samples=config.get("samples", "config/samples.tsv"),
        config="config/config.yaml",
        script="scripts/transmission.py",
    output:
        pairs="report/tables/snp_pairs.tsv",
        clusters="report/tables/clusters.tsv",
        heatmap="report/figures/snp_heatmap.png",
        timeline="report/figures/timeline.png",
        tree="report/figures/tree.png",
    log:
        "logs/transmission.log",
    conda:
        "../envs/plot.yaml"
    shell:
        "python {input.script} {input.matrix} {input.tree} {input.samples} {input.config} > {log} 2>&1"
