#!/bin/bash -l
#$ -cwd
#$ -l h_data=120G,h_rt=24:00:00,highp
#$ -j y
#$ -o ./job_out

export PATH=~/project-pasaniuc/software/anaconda3/bin:$PATH
export PYTHONNOUSERSITE=True

. /u/local/Modules/default/init/modules.sh

snakemake --cores 1
