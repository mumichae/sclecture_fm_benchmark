#!/usr/bin/env bash
set -e -x

export TMPDIR="/vol/data/tmp"
pipeline="$(realpath ../scAtlasTb)"

snakemake \
  --profile .profiles/local \
  --configfile \
    configs/subset/defaults_subset.yaml \
    configs/subset/qc_subset.yaml \
    configs/subset/integration_benchmark_subset.yaml \
  --snakefile $pipeline/workflow/Snakefile \
  "$@"
