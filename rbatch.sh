#!/bin/bash -l
#SBATCH --job-name=rbatch_1
#SBATCH --array=1-8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --partition brc
#SBATCH --mem-per-cpu=30000MB
#SBATCH --mail-type=END,FAIL
# script to submit an r batch job

number=$SLURM_ARRAY_TASK_ID
module load apps/R/3.6.0
Rscript 08_shearwaterML.R $number

