#!/usr/bin/env bash
set -e -x

pipeline="$(realpath ../sc-atlasing-toolbox)"

snakemake \
  --profile .profiles/local \
  --configfile \
    configs/example_config.yaml \
  --snakefile $pipeline/workflow/Snakefile \
  --use-conda \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds \
    $@
