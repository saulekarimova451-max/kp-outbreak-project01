# Project notes: what I did and why

Working log for the report and the defence. One section per phase, in plain words.
Numbers marked **[measure]** are filled in once the tool has run.

---

## Phase 0: Setup

- **Linux (WSL2 Ubuntu) on Windows.** Bioinformatics tools only run on Linux.
- **conda (Miniforge), channels conda-forge + bioconda, strict priority.** Installs exact tool versions in one command, so anyone can rebuild the same setup.
- **Six environments** (`kp`, `snippy`, `phylo`, `asm`, `amr`, `clair3`), recipes saved in `workflow/envs/`. Tools like Snippy, Gubbins and Clair3 need different library versions and break each other if installed together.
- **Git + GitHub with SSH.** Every phase is committed, so the history shows how the work progressed. Data files are excluded by `.gitignore`; scripts that download them are committed instead.
- **Computer:** 16 CPU threads, 7.6 GB RAM visible to Linux, about 950 GB free disk.

## Phase 1: The question and the threshold (Task 1)

- **Question:** were the *K. pneumoniae* isolates in this hospital spread from patient to patient (one strain, so an infection-control failure), or acquired separately?
- **Why genomes can answer it:** bacteria gain only a few mutations (SNPs) per year. Isolates from one chain of transmission over weeks differ by very few SNPs. Unrelated isolates of the same type differ by many more.
- **Threshold, fixed before any analysis:** ≤ 21 SNPs = genomically linked. Source: David et al. 2019, *Nature Microbiology*. They chose 21 because it best separated isolates from the same hospital from isolates from different hospitals, across 32 European countries. 22 to 50 SNPs = grey zone (needs epidemiological support); > 50 = unrelated.
- **Proof it was set in advance:** `config/config.yaml` was committed before any reads were downloaded (see `git log`).
- **Caveat:** the same resistance plasmid can move between different strains, so shared resistance genes alone do not prove transmission. Chromosome SNPs do.

## Phase 2: Data and ID reconciliation (Task 2)

### Dataset choice

- **First choice rejected: PRJNA551327** (Shanghai ICU, Chen et al. 2019). ENA showed 9 Illumina NovaSeq runs and **no Nanopore runs**, although the paper describes Nanopore sequencing. 4 of the 9 runs are labelled SINGLE but look like paired reads (about 100 bp per read, twice the read count). Without Nanopore data, Task 5 is impossible, so the dataset was dropped on day 1.
- **Chosen: PRJNA1251496.** Birmingham Heartlands Hospital (UK), carbapenem-resistant *K. pneumoniae* ST395, 2024. 8 isolates, each sequenced on **both** Illumina and Oxford Nanopore. Total download 12.15 GB.

### Retrieval (all scripted, all repeatable)

| Archive | Script / method | What it gave |
|---|---|---|
| ENA | `scripts/fetch_metadata.py` | run table: platform, model, read and base counts, file links, md5 |
| NCBI BioSample | `scripts/build_samples.py` | collection date, isolation source, host, location |
| NCBI Pathogen Detection | Isolates Browser export (`report/tables/pathogen_detection_isolates.tsv`) | SNP cluster, NCBI's AMR genotypes, assembly accessions |
| NCBI Assembly | `datasets summary genome` | assembly level and size per isolate |
| ENA (reads) | `scripts/fetch_reads.py` | FASTQ files, each verified by md5 |

### Problems found (full list in `report/tables/id_issues.tsv`)

- **Split BioSamples (8 of 8 isolates):** each isolate has two BioSample IDs, one for Illumina (SAMN48407xxx) and one for Nanopore (SAMN48004xxx). They could not be joined on BioSample.
- **Name suffixes (8 of 8):** Illumina samples end in `_trimmed`, Nanopore in `_final`. I stripped these labels and matched isolates by the shared name. All 8 matched.
- **Host is "restricted access" (8 of 8):** patient identity is hidden for privacy, so I **cannot tell which isolates come from the same patient**.
- **Two isolates share one sample number:** `0327436-KPC` (14 May) and `0327436-NDM` (8 May). Probably one patient carrying two strains. Kept as separate isolates.
- **"Create date" in Pathogen Detection is 2025** (upload date), not the 2024 collection date.
- **Lab labels vs. NCBI-computed genes disagree:** e.g. 6257373 labelled "NDM-4;NDM-5", NCBI computed blaNDM-1 + blaNDM-5.
- **The isolate labelled NDM has no carbapenemase gene in NCBI's results,** and its assembly is about 300 kb smaller (5.43 Mb vs 5.71 to 5.79 Mb). Hypothesis: it lost the plasmid carrying NDM. To test with Nanopore (Phase 6).
- **Download corruption:** one 2.3 GB file arrived at full size but failed the md5 check (gzip CRC error). The download script now deletes such files and re-downloads them (up to 3 attempts).

### What the metadata already shows

- 7 isolates collected 26 Apr to 14 May 2024 (under 3 weeks), 1 isolate on 10 Jul 2024. **Span about 11 weeks**, like the handbook scenario.
- 6 rectal swabs (screening, i.e. carriage) and 2 urine samples (likely infection).
- NCBI Pathogen Detection puts **all 8 in one SNP cluster (PDS000130824)**, together with 117 isolates from elsewhere. The closest pair differs by 1 SNP. NCBI's clustering is looser than my 21-SNP rule, so my analysis asks the sharper question: which pairs are close enough to mean direct transmission?
- All 8 NCBI assemblies are **contig level** (63 to 97 pieces, Illumina only). No complete genome exists, so I will build a complete reference from this outbreak's own Nanopore data.

## Phase 3: Platform reasoning (Task 3) — draft

| Property | Illumina | Oxford Nanopore |
|---|---|---|
| Instrument (ENA metadata) | HiSeq 4000, paired-end | GridION |
| Read length | about 241 bases per read on average (483 bases per pair ÷ 2), reads already trimmed by the submitters **[measure: seqkit stats]** | about 840 to 890 bases per read on average **[measure: NanoPlot N50]** |
| Coverage | 107x to 150x | 93x to 395x |
| Main error type | substitutions (one wrong letter); quality drops towards read ends | insertions and deletions, mostly in homopolymers (runs like AAAAA) |
| Good for | precise SNP calling | assembling whole chromosomes and plasmids, finding where resistance genes sit |
| Weak at | repeats and plasmids, so assemblies break into 60 to 100 pieces | single-base accuracy, especially in homopolymers |

**Metadata inconsistency to report:** ENA says HiSeq 4000, which produces reads up to 2 × 150 bp. But the reads average about 241 bases each, and the paper describes paired 250 bp reads. So the instrument label and the data disagree. The read lengths suggest a 2 × 250 bp run, which HiSeq 4000 cannot produce.

**Surprise to check:** a mean Nanopore read length of about 0.85 kb is short for long-read sequencing (often 5 to 20 kb). If NanoPlot confirms it, plasmid assembly may be harder than expected.

**48-hour turnaround (show the arithmetic in the report):**
- Data needed: 5.5 Mb genome × 8 isolates × 100x ≈ 4.4 Gb.
- Nanopore: one GridION flow cell with 8 barcoded isolates can produce this within hours, and data streams in real time.
- Illumina: library preparation plus a run of about 1 to 3 days, often batched with other customers' samples at a sequencing provider.
- Plus my measured pipeline runtime **[measure: Phase 8 benchmarks]**.
- Conclusion: Nanopore can meet 48 hours; outsourced Illumina usually cannot. Illumina remains the more accurate SNP ruler, which Phase 6 tests directly.
