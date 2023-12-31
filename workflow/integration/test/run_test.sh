#!/usr/bin/env bash
set -e -x

WORKDIR=$(dirname $(dirname $0))
cd $WORKDIR

#--snakefile $WORKDIR/Snakefile
snakemake --configfile test/config.yaml --use-conda $@

conda run --live-stream -n scanpy python test/run_assertions.py
#conda run --live-stream -n scvi-tools python test/run_scvi-tools_model_loading.py
#conda run --live-stream -n scarches python test/run_scarches_model_loading.py
