#!/usr/bin/env bash
set -e -x

WORKDIR=$(dirname $(dirname $0))
cd $WORKDIR

#--snakefile $WORKDIR/Snakefile
snakemake --configfile test/config_no_gpu.yaml --use-conda $@

conda run -p $CONDA_PREFIX/../scanpy python test/run_assertions.py --live-stream