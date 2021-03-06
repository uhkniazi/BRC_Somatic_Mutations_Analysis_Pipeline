# File: 05_bwa_array_job.R
# Auth: umar.niazi@kcl.as.uk
# DESC: create a parameter file and shell script to run array job on hpc
# Date: 18/10/2021


## set variables and source libraries
source('header.R')

setwd(gcRemoteDir)
csFiles = list.files('BRC_Somatic_Mutations_Analysis_Pipeline/dataExternal/remote/Trimmed/', 
                     pattern = '*.fastq')

# set header variables 
cvShell = '#!/bin/bash -l'
cvJobName = '#SBATCH --job-name=bwa-array'
cvNodes = '#SBATCH --nodes=1'
cvProcessors = '#SBATCH --ntasks=4'
# format d-hh:mm:ss
cvRuntime = '#SBATCH --time=6-20:005:00'
cvPartition = '#SBATCH --partition brc'
# How much memory you need.
# --mem will define memory per node and
# --mem-per-cpu will define memory per CPU/core. Choose one of those.
cvMemoryReserve = '#SBATCH --mem-per-cpu=6000MB'
# Turn on mail notification. There are many possible self-explaining values:
# NONE, BEGIN, END, FAIL, ALL (including all aforementioned)
# For more values, check "man sbatch"
cvMail = '#SBATCH --mail-type=END,FAIL'
# set array job loop
length(csFiles)
cvArrayJob = '#SBATCH --array=1-8'

# set the directory names
cvInput = 'input/'
cvOutput = 'output/Aligned/'
cvGeneIndex = '/users/k1625253/scratch/old-scratch_rosalind-legacy-import_2020-01-28/Data/MetaData/GenomeIndex/hg38_bwa/hg38.fa'

setwd(gcswd)
# create a parameter file and shell script
dir.create('AutoScripts')
oFile.param = file('AutoScripts/bwa_param.txt', 'wt')

# get unique file name after removing forward / reverse marks
csFiles = unique(gsub('_\\d\\.fastq.gz', '', csFiles))


temp = sapply(seq_along(csFiles), function(x){
  lf = list(paste0(csFiles[x], '_1.fastq.gz'),
            paste0(csFiles[x], '_2.fastq.gz'))
  # write bwa parameters
  in.r1 = paste0(cvInput, lf[[1]])
  in.r2 = paste0(cvInput, lf[[2]])
  out.sam = paste0(cvOutput, lf[[1]], '.sam')
  p1 = paste(in.r1, in.r2, out.sam, sep=' ')
  writeLines(p1, oFile.param)
})

close(oFile.param)

oFile = file('AutoScripts/bwa-mem.sh', 'wt')

writeLines(c(cvShell, cvJobName, cvArrayJob, cvNodes, cvProcessors, 
             cvRuntime, cvPartition, cvMemoryReserve, cvMail), oFile)

writeLines(c('# Autogenerated script from bwa_array_job.R', paste('# date', date())), oFile)
writeLines(c('# make sure directory paths exist before running script'), oFile)
writeLines('\n\n', oFile)

# module load
writeLines(c('module load apps/bwa/0.7.17-singularity'), oFile)
writeLines('\n\n', oFile)

## write array job lines
writeLines("# Parse parameter file to get variables.
number=$SLURM_ARRAY_TASK_ID
paramfile=bwa_param.txt
 
inr1=`sed -n ${number}p $paramfile | awk '{print $1}'`
inr2=`sed -n ${number}p $paramfile | awk '{print $2}'`
outsam=`sed -n ${number}p $paramfile | awk '{print $3}'`

# 9. Run the program.", oFile)
p1 = paste('bwa mem -t 4', cvGeneIndex, '$inr1', '$inr2', '>', '$outsam', sep=' ')
com = paste(p1)

writeLines(com, oFile)
writeLines('\n\n', oFile)

close(oFile)
