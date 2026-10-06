# Phase 7: resistance genes with BLAST against CARD (Task 7).
# Queries: NCBI's Illumina assembly of every isolate + our Flye (Nanopore) assemblies.

AMR_SAMPLES = ISOLATES + [f"flye-{i}" for i in config["ont"]["assemble"]]


def assembly_accession(wc):
    """NCBI assembly accession of an isolate, found via its Illumina BioSample in the Pathogen Detection export."""
    pdt = pd.read_csv(config["pathogen_detection"], sep="\t", dtype=str)
    biosample = SAMPLES.loc[wc.isolate, "illumina_biosample"]
    return pdt.loc[pdt["BioSample"] == biosample, "Assembly"].iloc[0]


def amr_query(wc):
    if wc.sample.startswith("flye-"):
        return f"results/ont/flye/{wc.sample[len('flye-'):]}/assembly.fasta"
    return f"resources/assemblies/{wc.sample}.fasta"


rule fetch_assembly:
    """Download NCBI's existing Illumina assembly for an isolate (saves re-assembling 8 genomes)."""
    output:
        "resources/assemblies/{isolate}.fasta",
    params:
        acc=assembly_accession,
        zip="resources/assemblies/{isolate}.zip",
    log:
        "logs/fetch_assembly/{isolate}.log",
    conda:
        "../envs/ref.yaml"
    shell:
        "datasets download genome accession {params.acc} --include genome --filename {params.zip} > {log} 2>&1 && "
        "unzip -p {params.zip} '*.fna' > {output} && rm {params.zip}"


rule card_db:
    """Download the current CARD release and build a BLAST nucleotide database."""
    output:
        fasta="resources/card/nucleotide_fasta_protein_homolog_model.fasta",
        index="resources/card/aro_index.tsv",
        version="resources/card/VERSION",
        done="resources/card/card_nt.done",
    log:
        "logs/card_db.log",
    conda:
        "../envs/amr.yaml"
    shell:
        "mkdir -p resources/card && "
        "wget -q -O resources/card/card.tar.bz2 https://card.mcmaster.ca/latest/data && "
        "tar -xjf resources/card/card.tar.bz2 -C resources/card && "
        "(grep -o '\"_version\": *\"[^\"]*\"' resources/card/card.json | head -1 > {output.version}; true) && "
        "sed -i '/^>/{{s/|/__/g; s/ /_/g}}' {output.fasta} && "
        "makeblastdb -in {output.fasta} -dbtype nucl -out resources/card/card_nt > {log} 2>&1 && "
        "touch {output.done}"


rule blast_card:
    """BLAST every contig against CARD's protein-homolog models (nucleotide sequences)."""
    input:
        query=amr_query,
        db="resources/card/card_nt.done",
    output:
        "results/amr/blast/{sample}.tsv",
    log:
        "logs/blast_card/{sample}.log",
    benchmark:
        "benchmarks/blast_card/{sample}.tsv"
    threads: THREADS
    conda:
        "../envs/amr.yaml"
    shell:
        "blastn -query {input.query} -db resources/card/card_nt -evalue 1e-10 -perc_identity 80 "
        "-max_target_seqs 50 -num_threads {threads} "
        "-outfmt '6 qseqid sseqid pident length qstart qend sstart send slen evalue bitscore' "
        "> {output} 2> {log}"


rule amr_summary:
    """Filter hits, predict the antibiogram, compare with NCBI, locate carbapenemases on Flye contigs."""
    input:
        blast=expand("results/amr/blast/{sample}.tsv", sample=AMR_SAMPLES),
        aro="resources/card/aro_index.tsv",
        samples=config.get("samples", "config/samples.tsv"),
        ncbi=config["pathogen_detection"],
        script="scripts/amr_summary.py",
    output:
        genes="report/tables/amr_genes.tsv",
        antibiogram="report/tables/antibiogram_predicted.tsv",
        ncbi="report/tables/amr_vs_ncbi.tsv",
        figure="report/figures/amr_genes.png",
    log:
        "logs/amr_summary.log",
    conda:
        "../envs/plot.yaml"
    shell:
        "python {input.script} --blast results/amr/blast --aro {input.aro} --samples {input.samples} "
        "--ncbi {input.ncbi} --flye results/ont/flye --outdir report > {log} 2>&1"
