#!/usr/bin/env bash
# Lock the exact package versions of every environment Snakemake built, into workflow/envs/locked/.
# Run from the project folder after a full workflow run:   bash scripts/export_envs.sh
set -euo pipefail
mkdir -p workflow/envs/locked
# Snakemake keeps a copy of each environment file next to the environment: .snakemake/conda/<hash>.yaml
for y in $(ls -tr .snakemake/conda/*.yaml); do
  name=$(awk '/^name:/ {print $2}' "$y")
  prefix="${y%.yaml}_"
  [ -d "$prefix" ] || continue
  conda env export -p "$prefix" --no-builds | grep -v '^prefix:' | sed "s/^name: .*/name: ${name}/" \
    > "workflow/envs/locked/${name}.lock.yaml"
  echo "locked ${name}  ($(grep -c '^  - ' workflow/envs/locked/${name}.lock.yaml) packages)"
done
