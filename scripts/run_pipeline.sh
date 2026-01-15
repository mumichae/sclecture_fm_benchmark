#!/usr/bin/env bash
set -e -x

export TMPDIR="/vol/data/tmp"
pipeline="$(realpath ../scAtlasTb)"

snakemake \
  --profile .profiles/local \
  --configfile \
    configs/defaults.yaml  \
    configs/qc.yaml  \
  --snakefile $pipeline/workflow/Snakefile \
    $@
