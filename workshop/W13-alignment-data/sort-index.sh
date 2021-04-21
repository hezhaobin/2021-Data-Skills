#!/bin/bash
########################
# author: Bin He
# date: 2020-04-20
# title: sort and index the alignment file
# use: qsub sort-index.sh
########################
# ----------------------
#  Scheduler parameters
# ----------------------
#  Specify the queue to use
#$ -q BIO-INSTR
#  give a name to the job, default to script name
#$ -N sort-index
#  execute the job from the current working directory
#$ -cwd
#  Request 10 slots in a shared memory environment
#$ -pe smp 10 
#  Specify the directory to save the job logs
#$ -o log/
#  Send e-mail at end/abort of job
#$ -m ea
# ---------------- edit the following part ---------------
#  E-mail address to send to (change <HawkID> to yours)
#$ -M HawkID@uiowa.edu
#---------------------------------------------------------

# these are useful flags to set to make the code more robust to failure
# copied from Vince Buffalo's Bioinformatic Data Analysis book
set -eo pipefail

# load samtools
module load stack/2020.2
module load samtools

# convert the SAM file to BAM
echo "-------- SAM -> BAM -------"
samtools view -b -@ 10 SRR6900282.aln.sam -o SRR6900282.bam # -@ 10 tells samtools to use 10 threads

# sort the BAM file
echo "-------- Sorting BAM -------"
samtools sort -@ 10 -T /localscratch/aln.sorted SRR6900282.sorted.bam

# index the sorted BAM file
echo "-------- Indexing BAM -------"
samtools index -@ 10 SRR6900282.sorted.bam
