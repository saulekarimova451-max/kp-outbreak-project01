# *Klebsiella pneumoniae* hospital outbreak: an Illumina + Nanopore genomic investigation

Course project 01 ("The Outbreak Investigation"), Introduction to Bioinformatics, Astana IT University.
Author: *[your name]* (individual project, approved by the instructor).

## The question and the answer

**Were the carbapenem-resistant *K. pneumoniae* isolates from one hospital spread between patients, or acquired separately?**

| Finding | Evidence |
|---|---|
| **One outbreak.** All 8 isolates are genomically linked. | Every pair is 1 to 9 core-genome SNPs apart (median 5). Threshold, fixed before analysis: ≤ 21 SNPs (David et al. 2019). 28 of 28 pairs linked, 1 cluster. |
| **Two sub-lineages, still circulating in July.** | Tree: group A (earliest isolate 0326576 + patient 0327436, bootstrap 97) and a star-shaped group B around 0327127. The July isolate 6257373 is 3 SNPs from 0327127. |
| **NDM carbapenemase in 7 of 8 isolates**, usually two alleles (NDM-1 + NDM-5). | BLAST against CARD; agrees with NCBI AMRFinderPlus at gene-family level. |
| **One isolate lost its resistance plasmid.** | 0327436-NDM carries no carbapenemase and lacks sul1, dfrA12, mph(A), qnrB and aadA2; its genome is about 300 kb smaller. |
| **Nanopore alone would have missed the outbreak** with this data. | Nanopore-only SNP distances: 36 to 245; 0 of 28 pairs linked. The reads are short (N50 about 0.8 kb) and Q13.7. Illumina remains necessary for SNP-level transmission calls. |

Full interpretation: `report/` (report PDF, figures, tables). Running notes: `report/notes.md`.

## Data

All data are public and downloaded by scripts.

| What | Accession / source |
|---|---|
| Reads (Illumina + Nanopore, 8 isolates) | BioProject [PRJNA1251496](https://www.ebi.ac.uk/ena/browser/view/PRJNA1251496), runs SRR33529050–57 (Illumina), SRR33186507–14 (Nanopore) |
| Sample metadata | NCBI BioSample (two per isolate: SAMN48407216–23 Illumina, SAMN48004xxx Nanopore) |
| Outbreak cluster and AMR baseline | NCBI Pathogen Detection, SNP cluster PDS000130824 (`report/tables/pathogen_detection_isolates.tsv`, exported from the Isolates Browser) |
| Reference genome | GCA_050993895.1 (NCBI assembly of the earliest isolate, 0326576) |
| Illumina assemblies for AMR screening | GCA_050993855.1 … GCA_050994035.1 (one per isolate) |
| Resistance gene database | CARD, latest release at run time (version written to `resources/card/VERSION`) |

Dataset choice: PRJNA551327 (Shanghai ICU) was evaluated first and rejected because the Nanopore runs described in its paper were never deposited.

## Requirements

- Linux (tested on Ubuntu under WSL2 on Windows 11), 8 GB RAM minimum, about 30 GB free disk
- [Miniforge](https://github.com/conda-forge/miniforge) (conda), with channels `conda-forge` and `bioconda`, strict priority
- Internet access (ENA, NCBI, CARD)

```bash
conda create -n kp -c conda-forge -c bioconda snakemake pandas
conda activate kp
```

Every analysis tool is installed automatically by Snakemake into its own environment (`workflow/envs/*.yaml`; exact versions in `workflow/envs/locked/`).

## Quick test (about 20 to 40 minutes, mostly installing tools)

From a fresh clone, on the small dataset in `test/` (2 isolates, ~5x coverage):

```bash
git clone git@github.com:saulekarimova451-max/kp-outbreak-project01.git
cd kp-outbreak-project01
snakemake --sdm conda --cores 2 test --configfile config/test.yaml
```

Success = `results/qc/multiqc_report.html`, and Illumina and Nanopore variant files for both test isolates.

## Full analysis (one command after downloading the data)

```bash
# 1. metadata and reads (12 GB; resumable, md5-verified)
python scripts/fetch_metadata.py PRJNA1251496
python scripts/build_samples.py PRJNA1251496
python scripts/fetch_reads.py data/meta/PRJNA1251496_ena_runs.tsv data/raw

# 2. everything else: QC, SNPs, tree, Nanopore, resistance genes, validation, figures
snakemake --sdm conda --cores 4
```

Runtime on a 16-thread laptop with 4 cores in use: see `report/tables/runtime.tsv`. Single steps can be run by name:

| Target | What it does |
|---|---|
| `qc_illumina`, `qc_ont` | FastQC, fastp, MultiQC, seqkit; NanoPlot |
| `snps` | Snippy (BWA-MEM + freebayes) per isolate, core alignment, Gubbins, snp-dists, IQ-TREE |
| `figures` | transmission calls at the pre-registered threshold, heatmap, timeline, tree |
| `ont` | Filtlong, minimap2, Clair3, Illumina vs Nanopore comparison, Flye assemblies |
| `amr` | BLAST against CARD, predicted antibiogram, comparison with NCBI |
| `validate` | 100 known SNPs simulated into the reference: measured precision and recall of the SNP pipeline |
| `runtime` | runtime per step and a scale estimate for 100 and 1,000 isolates |

All thresholds are in `config/config.yaml`. The SNP thresholds were committed before any data were downloaded (`git log config/config.yaml`).

## Main outputs

| File | Content |
|---|---|
| `report/tables/snp_pairs.tsv`, `clusters.tsv` | every isolate pair: SNPs, days apart, linked / grey zone / unrelated |
| `report/figures/snp_heatmap.png`, `timeline.png`, `tree.png` | the outbreak picture |
| `report/tables/illumina_vs_nanopore.tsv`, `snp_matrix_nanopore.tsv`, `nanopore_only_snp_profile.tsv` | platform comparison |
| `report/tables/amr_genes.tsv`, `antibiogram_predicted.tsv`, `amr_vs_ncbi.tsv` | resistance genes and predicted antibiogram |
| `report/tables/validation_snippy.tsv` | measured SNP-calling accuracy |
| `report/tables/id_issues.tsv` | every identifier and metadata problem found, and how it was resolved |
| `report/tables/runtime.tsv` | measured runtime per step |

## Repository structure

```text
config/      config.yaml (thresholds), samples.tsv (isolate ↔ runs ↔ BioSamples), test.yaml
workflow/    Snakefile, rules/*.smk, envs/*.yaml (+ envs/locked/)
scripts/     data retrieval, analysis and figure scripts (Python, bash)
test/        tiny test dataset
report/      report, figures, tables, notes.md
benchmarks/  runtime of every job (written by Snakemake)
data/, results/, resources/   created by the workflow, not committed
```

## Reproducibility notes

- Environments: `workflow/envs/*.yaml` define each tool set; `scripts/export_envs.sh` writes the exact versions used to `workflow/envs/locked/`.
- Container: `Dockerfile` generated with `snakemake --containerize`; not built locally (no Docker on the development laptop).
- Benchmarks record runtime; memory is not recorded under WSL2 (Snakemake cannot read it there).
- One manual step: the NCBI Pathogen Detection table was exported from the web browser (no programmatic export is offered); the file is committed.

## Tools

fastp, FastQC, MultiQC, seqkit, NanoPlot, Filtlong, BWA-MEM, minimap2, samtools, mosdepth, Snippy, freebayes, Gubbins, snp-sites, snp-dists, IQ-TREE, Clair3, Flye, BLAST+, CARD, NCBI Datasets, Snakemake. Versions: `workflow/envs/locked/`.

AI assistance is disclosed in the report appendix.
