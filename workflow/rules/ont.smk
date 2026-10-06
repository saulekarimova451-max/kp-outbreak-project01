# Phase 6: Nanopore (TGS) analysis and comparison with Illumina (Task 5).
# Settings in config/config.yaml, section "ont" (chosen from the NanoPlot results).


rule faidx_reference:
    """Index the reference so tools can jump to any position (needed by Clair3)."""
    input:
        "resources/reference.fasta",
    output:
        "resources/reference.fasta.fai",
    conda:
        "../envs/ont.yaml"
    shell:
        "samtools faidx {input}"


rule filtlong:
    """Drop short and low-quality reads, then keep the best reads up to ~100x coverage."""
    input:
        ont_reads,
    output:
        "results/ont/filtered/{isolate}.fastq.gz",
    params:
        min_len=config["ont"]["min_length"],
        min_q=config["ont"]["min_mean_q"],
        target=config["ont"]["target_bases"],
    log:
        "logs/filtlong/{isolate}.log",
    benchmark:
        "benchmarks/filtlong/{isolate}.tsv"
    conda:
        "../envs/ont.yaml"
    shell:
        "filtlong --min_length {params.min_len} --min_mean_q {params.min_q} "
        "--target_bases {params.target} {input} 2> {log} | gzip > {output}"


rule map_ont:
    """Map Nanopore reads with minimap2 (built for long, error-prone reads), sort and index."""
    input:
        reads="results/ont/filtered/{isolate}.fastq.gz",
        ref="resources/reference.fasta",
    output:
        bam="results/ont/bam/{isolate}.bam",
        bai="results/ont/bam/{isolate}.bam.bai",
    log:
        "logs/map_ont/{isolate}.log",
    benchmark:
        "benchmarks/map_ont/{isolate}.tsv"
    threads: THREADS
    conda:
        "../envs/ont.yaml"
    shell:
        "(minimap2 -ax map-ont -t {threads} {input.ref} {input.reads} "
        "| samtools sort -@ 2 -o {output.bam} - && samtools index {output.bam}) 2> {log}"


rule clair3:
    """Haploid variant calling with a neural network trained on bacterial R10.4.1 Nanopore data."""
    input:
        bam="results/ont/bam/{isolate}.bam",
        ref="resources/reference.fasta",
        fai="resources/reference.fasta.fai",
    output:
        "results/ont/clair3/{isolate}/merge_output.vcf.gz",
    params:
        outdir="results/ont/clair3/{isolate}",
        model=config["ont"]["clair3_model"],
    log:
        "logs/clair3/{isolate}.log",
    benchmark:
        "benchmarks/clair3/{isolate}.tsv"
    threads: THREADS
    conda:
        "../envs/clair3.yaml"
    shell:
        "run_clair3.sh --bam_fn={input.bam} --ref_fn={input.ref} --threads={threads} "
        "--platform=ont --model_path=$CONDA_PREFIX/bin/models/{params.model} "
        "--output={params.outdir} --haploid_precise --include_all_ctgs --no_phasing_for_fa "
        "> {log} 2>&1"


rule compare_platforms:
    """Illumina (Snippy) vs Nanopore (Clair3) variant calls for every isolate with both."""
    input:
        illumina=expand("results/snippy/{isolate}/snps.vcf", isolate=ISOLATES),
        ont=expand("results/ont/clair3/{isolate}/merge_output.vcf.gz", isolate=ISOLATES),
        ref="resources/reference.fasta",
        script="scripts/compare_platforms.py",
    output:
        summary="report/tables/illumina_vs_nanopore.tsv",
        discordant="report/tables/illumina_vs_nanopore_discordant.tsv",
        figure="report/figures/illumina_vs_nanopore.png",
    params:
        isolates=",".join(ISOLATES),
    log:
        "logs/compare_platforms.log",
    conda:
        "../envs/plot.yaml"
    shell:
        "python {input.script} --ref {input.ref} --isolates {params.isolates} "
        "--illumina 'results/snippy/{{iso}}/snps.vcf' "
        "--ont 'results/ont/clair3/{{iso}}/merge_output.vcf.gz' "
        "--summary {output.summary} --discordant {output.discordant} --figure {output.figure} "
        "> {log} 2>&1"


rule flye:
    """Long-read assembly (only for isolates listed in config ont.assemble)."""
    input:
        "results/ont/filtered/{isolate}.fastq.gz",
    output:
        fasta="results/ont/flye/{isolate}/assembly.fasta",
        info="results/ont/flye/{isolate}/assembly_info.txt",
    params:
        outdir="results/ont/flye/{isolate}",
    log:
        "logs/flye/{isolate}.log",
    benchmark:
        "benchmarks/flye/{isolate}.tsv"
    threads: THREADS
    conda:
        "../envs/asm.yaml"
    shell:
        "flye --nano-hq {input} --out-dir {params.outdir} --threads {threads} > {log} 2>&1"


rule nanopore_distances:
    """Would Nanopore alone give the same outbreak call? Distances + profile of Nanopore-only SNPs."""
    input:
        illumina=expand("results/snippy/{isolate}/snps.vcf", isolate=ISOLATES),
        ont=expand("results/ont/clair3/{isolate}/merge_output.vcf.gz", isolate=ISOLATES),
        ref="resources/reference.fasta",
        script="scripts/nanopore_distances.py",
        helper="scripts/compare_platforms.py",
    output:
        matrix="report/tables/snp_matrix_nanopore.tsv",
        profile="report/tables/nanopore_only_snp_profile.tsv",
        figure="report/figures/illumina_vs_nanopore_distances.png",
    params:
        isolates=",".join(ISOLATES),
    log:
        "logs/nanopore_distances.log",
    conda:
        "../envs/plot.yaml"
    shell:
        "python {input.script} --ref {input.ref} --isolates {params.isolates} "
        "--illumina 'results/snippy/{{iso}}/snps.vcf' "
        "--ont 'results/ont/clair3/{{iso}}/merge_output.vcf.gz' "
        "--matrix {output.matrix} --profile {output.profile} --figure {output.figure} > {log} 2>&1"
