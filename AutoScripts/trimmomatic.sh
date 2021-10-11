#!/bin/bash -l
#SBATCH --job-name=trim-array
#SBATCH --array=1-8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-05:005:00
#SBATCH --partition brc
#SBATCH --mem-per-cpu=10000MB
#SBATCH --mail-type=END,FAIL
# Autogenerated script from write_trimmomatic_script.R
# date Mon Oct 11 16:44:12 2021
# make sure directory paths exist before running script



module load apps/trimmomatic/0.39



# Parse parameter file to get variables.
number=$SLURM_ARRAY_TASK_ID
paramfile=trimmomatic_param.txt
 
inr1=`sed -n ${number}p $paramfile | awk '{print $1}'`
inr2=`sed -n ${number}p $paramfile | awk '{print $2}'`
outr1=`sed -n ${number}p $paramfile | awk '{print $3}'`
outr1up=`sed -n ${number}p $paramfile | awk '{print $4}'`
outr2=`sed -n ${number}p $paramfile | awk '{print $5}'`
outr2up=`sed -n ${number}p $paramfile | awk '{print $6}'`

# 9. Run the program.
trimmomatic PE $inr1 $inr2 $outr1 $outr1up $outr2 $outr2up ILLUMINACLIP:/users/k1625253/scratch/old-scratch_rosalind-legacy-import_2020-01-28/Data/MetaData/trimmomatic_adapters.fa:2:30:10:8:true



