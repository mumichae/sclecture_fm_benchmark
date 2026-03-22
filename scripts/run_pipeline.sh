#!/usr/bin/env bash
set -e -x

export TMPDIR="/vol/data/tmp"
pipeline="$(realpath ../scAtlasTb)"

snakemake \
  --profile .profiles/local \
  --configfile \
    configs/defaults.yaml  \
    configs/marker_genes.yaml  \
    configs/integration.yaml  \
  --snakefile $pipeline/workflow/Snakefile \
    $@
