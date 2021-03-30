#!/bin/bash
#./bash/RNAseq-setup-start.sh
#purpose: copy read files from BOX into Pegasus scratch space and rename for processing
#To start this job from the anti_phys directory, use:
#bsub -P transcriptomics < ./bash/wget_box_reads.sh

#BSUB -J wget_box_reads
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o wget_box_reads%J.out
#BSUB -e wget_box_reads%J.err
#BSUB -n 8
#BSUB -u mconnelly@rsmas.miami.edu
#BSUB -N

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/projects/transcriptomics/mikeconnelly"
prodir="/scratch/projects/transcriptomics/mikeconnelly/projects/coral_antibiotics_physiology"

wget -L https://miami.box.com/s/oc47nri78qo1snu4hvahqk3sjhnjmapq/*.gz 
