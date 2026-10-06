# Phase 4: quality control (Task 4) and read-length facts for Task 3.
# Thresholds come from config/config.yaml (section "qc"), fixed before analysis.


rule fastqc_raw:
    """FastQC on the reads exactly as downloaded."""
    input:
        unpack(illumina_reads),
    output:
        directory("results/qc/fastqc_raw/{isolate}"),
    log:
        "logs/fastqc_raw/{isolate}.log",
    benchmark:
        "benchmarks/fastqc_raw/{isolate}.tsv"
    threads: 2
    conda:
        "../envs/qc.yaml"
    shell:
        "mkdir -p {output} && fastqc -t {threads} --outdir {output} {input.r1} {input.r2} > {log} 2>&1"


rule fastp:
    """Trim adapters and low-quality 3' ends, drop short reads, remove duplicates."""
    input:
        unpack(illumina_reads),
    output:
        r1="results/trimmed/{isolate}_R1.fastq.gz",
        r2="results/trimmed/{isolate}_R2.fastq.gz",
        json="results/qc/fastp/{isolate}.fastp.json",
        html="results/qc/fastp/{isolate}.fastp.html",
    params:
        q=config["qc"]["min_base_quality"],
        minlen=config["qc"]["min_read_length_illumina"],
    log:
        "logs/fastp/{isolate}.log",
    benchmark:
        "benchmarks/fastp/{isolate}.tsv"
    threads: 4
    conda:
        "../envs/qc.yaml"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "--detect_adapter_for_pe "
        "--cut_right --cut_right_window_size 4 --cut_right_mean_quality {params.q} "
        "--length_required {params.minlen} "
        "--dedup "
        "--thread {threads} "
        "--report_title {wildcards.isolate} "
        "--json {output.json} --html {output.html} > {log} 2>&1"


rule fastqc_trimmed:
    """FastQC again after fastp, to show what trimming changed."""
    input:
        r1="results/trimmed/{isolate}_R1.fastq.gz",
        r2="results/trimmed/{isolate}_R2.fastq.gz",
    output:
        directory("results/qc/fastqc_trimmed/{isolate}"),
    log:
        "logs/fastqc_trimmed/{isolate}.log",
    benchmark:
        "benchmarks/fastqc_trimmed/{isolate}.tsv"
    threads: 2
    conda:
        "../envs/qc.yaml"
    shell:
        "mkdir -p {output} && fastqc -t {threads} --outdir {output} {input.r1} {input.r2} > {log} 2>&1"


rule seqkit_illumina:
    """Read counts and length statistics, raw vs trimmed (measured read length for Task 3)."""
    input:
        raw=RAW_ILLUMINA,
        trimmed=expand("results/trimmed/{isolate}_R{m}.fastq.gz", isolate=ISOLATES, m=[1, 2]),
    output:
        "results/qc/seqkit_illumina.tsv",
    log:
        "logs/seqkit_illumina.log",
    benchmark:
        "benchmarks/seqkit_illumina.tsv"
    threads: THREADS
    conda:
        "../envs/qc.yaml"
    shell:
        "seqkit stats --all --tabular --threads {threads} {input.raw} {input.trimmed} > {output} 2> {log}"


rule multiqc:
    """One combined report of FastQC (before and after) and fastp."""
    input:
        expand("results/qc/fastqc_raw/{isolate}", isolate=ISOLATES),
        expand("results/qc/fastqc_trimmed/{isolate}", isolate=ISOLATES),
        expand("results/qc/fastp/{isolate}.fastp.json", isolate=ISOLATES),
    output:
        "results/qc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    conda:
        "../envs/qc.yaml"
    shell:
        "multiqc --force --outdir results/qc --filename multiqc_report.html "
        "results/qc/fastqc_raw results/qc/fastqc_trimmed results/qc/fastp > {log} 2>&1"


rule nanoplot:
    """Nanopore read length and quality distribution (raw reads)."""
    input:
        ont_reads,
    output:
        "results/qc/nanoplot/{isolate}/NanoStats.txt",
    params:
        outdir="results/qc/nanoplot/{isolate}",
    log:
        "logs/nanoplot/{isolate}.log",
    benchmark:
        "benchmarks/nanoplot/{isolate}.tsv"
    threads: 4
    conda:
        "../envs/qc.yaml"
    shell:
        "NanoPlot --fastq {input} --outdir {params.outdir} --threads {threads} "
        "--tsv_stats --loglength > {log} 2>&1"
