#!/usr/bin/env bash
set -e -x

snakemake exploration_all integration_all \
  --profile .profiles/czbiohub \
  --configfile \
    configs/computational_resources/czbiohub.yaml  \
    configs/integration/config.yaml \
    $@
