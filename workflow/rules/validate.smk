# Phase 5.3 + 8: measured accuracy of the Illumina SNP pipeline, and runtime summary.

VALIDATION = config.get("validation", {"n_snps": 100, "seed": 42, "coverage": 100, "read_length": 241})


rule mutate_reference:
    """Reference copy with N known random SNPs (the ground truth)."""
    input:
        ref="resources/reference.fasta",
        script="scripts/mutate_reference.py",
    output:
        fasta="results/validation/mutated.fasta",
        truth="results/validation/truth.tsv",
    params:
        n=VALIDATION["n_snps"],
        seed=VALIDATION["seed"],
    conda:
        "../envs/plot.yaml"
    shell:
        "python {input.script} {input.ref} {output.fasta} {output.truth} {params.n} {params.seed}"


rule simulate_reads:
    """Simulated Illumina pairs from the mutated genome (wgsim: 0.1% base errors, no extra mutations)."""
    input:
        "results/validation/mutated.fasta",
    output:
        r1="results/validation/sim_R1.fastq.gz",
        r2="results/validation/sim_R2.fastq.gz",
    params:
        pairs=lambda wc: int(VALIDATION["coverage"] * 5_714_824 / (2 * VALIDATION["read_length"])),
        rl=VALIDATION["read_length"],
        seed=VALIDATION["seed"],
    log:
        "logs/simulate_reads.log",
    conda:
        "../envs/sim.yaml"
    shell:
        "wgsim -N {params.pairs} -1 {params.rl} -2 {params.rl} -e 0.001 -r 0 -R 0 -S {params.seed} "
        "{input} results/validation/sim_R1.fastq results/validation/sim_R2.fastq > {log} 2>&1 && "
        "gzip -f results/validation/sim_R1.fastq results/validation/sim_R2.fastq"


rule snippy_validation:
    """The same Snippy settings as the real isolates, run on the simulated reads."""
    input:
        r1="results/validation/sim_R1.fastq.gz",
        r2="results/validation/sim_R2.fastq.gz",
        ref="resources/reference.fasta",
    output:
        "results/validation/snippy/snps.vcf",
    params:
        mincov=config["snippy"]["mincov"],
        minfrac=config["snippy"]["minfrac"],
    log:
        "logs/snippy_validation.log",
    benchmark:
        "benchmarks/validation_snippy.tsv"
    threads: THREADS
    conda:
        "../envs/snippy.yaml"
    shell:
        "snippy --force --cpus {threads} --outdir results/validation/snippy --ref {input.ref} "
        "--R1 {input.r1} --R2 {input.r2} --mincov {params.mincov} --minfrac {params.minfrac} > {log} 2>&1"


rule score_validation:
    """Precision, recall and F1 of the SNP calls against the known truth."""
    input:
        truth="results/validation/truth.tsv",
        vcf="results/validation/snippy/snps.vcf",
        bench="benchmarks/validation_snippy.tsv",
        script="scripts/score_validation.py",
    output:
        "report/tables/validation_snippy.tsv",
    log:
        "logs/score_validation.log",
    conda:
        "../envs/plot.yaml"
    shell:
        "python {input.script} {input.truth} {input.vcf} {input.bench} {output} > {log} 2>&1"


rule runtime_summary:
    """Runtime per step from all benchmark files, plus a scale estimate (run last)."""
    input:
        SNPS + ONT + AMR + FIGURES + ["report/tables/validation_snippy.tsv"],
        script="scripts/runtime_summary.py",
    output:
        table="report/tables/runtime.tsv",
        figure="report/figures/runtime.png",
    params:
        n=len(ISOLATES),
        cores=THREADS,
    log:
        "logs/runtime_summary.log",
    conda:
        "../envs/plot.yaml"
    shell:
        "python {input.script} benchmarks {output.table} {output.figure} {params.n} {params.cores} > {log} 2>&1"
