#!/usr/bin/env bash
# Build the small test dataset in test/ (about 40 MB): 2 isolates, Illumina + Nanopore, ~5x coverage.
# Run once from the project folder, with the real reads already in data/raw/:
#     conda activate kp && bash scripts/make_test_data.sh
# The test then runs from a fresh clone with:  snakemake --sdm conda --cores 2 test --configfile config/test.yaml
set -euo pipefail
ISOLATES=("0326575" "6257373")
COV=5
mkdir -p test/data
head -1 config/samples.tsv > test/samples.tsv
col() { head -1 config/samples.tsv | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1; }
ILL_COL=$(col illumina_run); ONT_COL=$(col ont_run)
for iso in "${ISOLATES[@]}"; do
  row=$(awk -F'\t' -v i="$iso" '$1==i' config/samples.tsv)
  echo "$row" >> test/samples.tsv
  ill=$(echo "$row" | cut -f"$ILL_COL")
  ont=$(echo "$row" | cut -f"$ONT_COL")
  rasusa reads --coverage "$COV" --genome-size 5.5mb --seed 1 \
    -o "test/data/${ill}_1.fastq.gz" -o "test/data/${ill}_2.fastq.gz" \
    "data/raw/${ill}_1.fastq.gz" "data/raw/${ill}_2.fastq.gz"
  rasusa reads --coverage "$COV" --genome-size 5.5mb --seed 1 \
    -o "test/data/${ont}_1.fastq.gz" "data/raw/${ont}_1.fastq.gz"
done
ls -lh test/data
du -sh test
