#!/bin/bash
#./bash/tagseq/unique_featurecounts_symC1.sh
#purpose: mark duplicate reads and quantify STAR-aligned BAM files at the gene level using featureCounts
#To start this job from the anti_phys directory, use:
#bsub -P transcriptomics < ./bash/tagseq/unique_featurecounts_symC1.sh

#BSUB -J unique_featurecounts_symC1
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o unique_featurecounts_symC1%J.out
#BSUB -e unique_featurecounts_symC1%J.err
#BSUB -n 8
#BSUB -u mconnelly@rsmas.miami.edu
#BSUB -N

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/projects/transcriptomics/mikeconnelly"
prodir="/scratch/projects/transcriptomics/mikeconnelly/projects/anti_phys"

# making a list of sample names
samples=$(cat ${prodir}/data/quantseq_samples.txt)

module load samtools
module load subread

#for samp in $samples ; do
#samtools index ${prodir}/outputs/alignments/${samp}_Aligned.sortedByCoord.out.bam
#samtools view -q 255 -Sub ${prodir}/outputs/STARalign_Pdam/${samp}_PdamAligned.sortedByCoord.out.bam \
#-o ${prodir}/outputs/STARalign_Pdam/${samp}_Aligned.sortedByCoord.out.uniq.bam
#done

featureCounts -p -T 8 -t gene \
-g ID \
-a ${mcs}/sequences/genomes/symbiodinium/symC1_genome.gff \
-o ${prodir}/outputs/STARcounts_SymC1/PocAnti_SymC1.counts \
${prodir}/outputs/STARalign_Pdam/*uniq.bam
