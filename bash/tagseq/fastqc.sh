#!/bin/bash
#./bash/fastqc.sh
#purpose: quality checking of raw RNAseq reads using FASTQC on Pegasus compute node
#To start this job from the anti_phys directory, use:
#bsub -P transcriptomics < /scratch/projects/transcriptomics/mikeconnelly/projects/anti_phys/bash/fastqc.sh

#BSUB -J fastqc
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o fastqc%J.out
#BSUB -e fastqc%J.err
#BSUB -n 12
#BSUB -u mconnelly@rsmas.miami.edu
#BSUB -N

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/projects/transcriptomics/mikeconnelly"
prodir="/scratch/projects/transcriptomics/mikeconnelly/projects/anti_phys"

module load java/1.8.0_60
module load fastqc/0.10.1
fastqc \
${prodir}/outputs/trimmed_reads/* \
--outdir ${prodir}/outputs/trimqcs/
