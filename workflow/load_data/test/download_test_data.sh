#!/usr/bin/env bash
set -e -x

WORKDIR=$(dirname $(dirname $0))
cd $WORKDIR

#--snakefile $WORKDIR/Snakefile
snakemake test/out/load_data/download/SchulteSchrepping2020.h5ad --configfile test/config.yaml --use-conda --printshellcmds $@
cp test/out/load_data/download/SchulteSchrepping2020.h5ad test/data/SchulteSchrepping2020.h5ad
