# Phase 5: reference, mapping, variant calling, SNP distances, tree (Task 6).


rule fetch_reference:
    """Download the outbreak reference genome from NCBI (accession set in config.yaml)."""
    output:
        "resources/reference.fasta",
    params:
        acc=config["reference"]["accession"],
    log:
        "logs/fetch_reference.log",
    conda:
        "../envs/ref.yaml"
    shell:
        "datasets download genome accession {params.acc} --include genome "
        "--filename resources/reference.zip > {log} 2>&1 && "
        "unzip -p resources/reference.zip '*.fna' > {output} && rm resources/reference.zip"


rule snippy:
    """Map trimmed Illumina reads to the reference (BWA-MEM) and call haploid variants (freebayes)."""
    input:
        r1="results/trimmed/{isolate}_R1.fastq.gz",
        r2="results/trimmed/{isolate}_R2.fastq.gz",
        ref="resources/reference.fasta",
    output:
        vcf="results/snippy/{isolate}/snps.vcf",
        bam="results/snippy/{isolate}/snps.bam",
        aln="results/snippy/{isolate}/snps.aligned.fa",
    params:
        outdir="results/snippy/{isolate}",
        mincov=config["snippy"]["mincov"],
        minfrac=config["snippy"]["minfrac"],
    log:
        "logs/snippy/{isolate}.log",
    benchmark:
        "benchmarks/snippy/{isolate}.tsv"
    threads: THREADS
    conda:
        "../envs/snippy.yaml"
    shell:
        "snippy --force --cpus {threads} --outdir {params.outdir} --ref {input.ref} "
        "--R1 {input.r1} --R2 {input.r2} --mincov {params.mincov} --minfrac {params.minfrac} "
        "> {log} 2>&1"


rule coverage:
    """Mean depth and breadth of coverage per isolate (QC thresholds in config.yaml)."""
    input:
        bam="results/snippy/{isolate}/snps.bam",
    output:
        summary="results/coverage/{isolate}.mosdepth.summary.txt",
        thresholds="results/coverage/{isolate}.thresholds.bed.gz",
    params:
        prefix="results/coverage/{isolate}",
    log:
        "logs/coverage/{isolate}.log",
    threads: 4
    conda:
        "../envs/coverage.yaml"
    shell:
        "samtools index {input.bam} 2> {log} && "
        "mosdepth --threads {threads} --no-per-base --by 100000 --thresholds 1,10,30 "
        "{params.prefix} {input.bam} >> {log} 2>&1"


rule snippy_core:
    """Core-genome alignment: positions covered in every isolate."""
    input:
        dirs=expand("results/snippy/{isolate}/snps.vcf", isolate=ISOLATES),
        ref="resources/reference.fasta",
    output:
        full="results/core/core.full.aln",
        clean="results/core/clean.full.aln",
        snps="results/core/core.aln",
        stats="results/core/core.txt",
    params:
        dirs=expand("results/snippy/{isolate}", isolate=ISOLATES),
    log:
        "logs/snippy_core.log",
    conda:
        "../envs/snippy.yaml"
    shell:
        "snippy-core --ref {input.ref} --prefix results/core/core {params.dirs} > {log} 2>&1 && "
        "snippy-clean_full_aln {output.full} > {output.clean} 2>> {log}"


rule snp_matrix_raw:
    """Pairwise SNP distances from the core SNP alignment (before recombination filtering)."""
    input:
        "results/core/core.aln",
    output:
        "results/core/snp_matrix_raw.tsv",
    conda:
        "../envs/snippy.yaml"
    shell:
        "snp-dists -b {input} > {output}"


rule gubbins:
    """Find and mask recombined regions, which would inflate SNP distances."""
    input:
        "results/core/clean.full.aln",
    output:
        filtered="results/core/gubbins.filtered_polymorphic_sites.fasta",
    log:
        "logs/gubbins.log",
    benchmark:
        "benchmarks/gubbins.tsv"
    threads: THREADS
    conda:
        "../envs/phylo.yaml"
    shell:
        "cd results/core && run_gubbins.py --threads {threads} --prefix gubbins clean.full.aln "
        "> ../../{log} 2>&1"


rule snp_matrix_clean:
    """SNP distances after removing recombination: the matrix used for the transmission call."""
    input:
        "results/core/gubbins.filtered_polymorphic_sites.fasta",
    output:
        aln="results/core/clean.core.aln",
        matrix="results/core/snp_matrix_clean.tsv",
    conda:
        "../envs/phylo.yaml"
    shell:
        "snp-sites -c -o {output.aln} {input} && snp-dists -b {output.aln} > {output.matrix}"


rule iqtree:
    """Maximum-likelihood tree from the clean SNP alignment, with 1000 ultrafast bootstraps."""
    input:
        "results/core/clean.core.aln",
    output:
        "results/tree/outbreak.treefile",
    log:
        "logs/iqtree.log",
    benchmark:
        "benchmarks/iqtree.tsv"
    threads: 4
    conda:
        "../envs/phylo.yaml"
    shell:
        "$(command -v iqtree3 || command -v iqtree2 || command -v iqtree) -s {input} -m GTR+ASC -B 1000 -T {threads} --prefix results/tree/outbreak -redo "
        "> {log} 2>&1"
