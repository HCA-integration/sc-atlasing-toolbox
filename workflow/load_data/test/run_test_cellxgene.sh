#!/usr/bin/env bash
set -e -x

WORKDIR=$(dirname $(dirname $0))
cd $WORKDIR

#--snakefile $WORKDIR/Snakefile
snakemake --configfile test/configs/cellxgene.yaml --rerun-incomplete --use-conda --printshellcmds $@

conda run -n scanpy python test/run_assertions.py --live-stream