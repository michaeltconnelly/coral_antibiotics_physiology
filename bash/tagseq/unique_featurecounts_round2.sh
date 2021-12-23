#!/bin/bash
#./bash/tagseq/unique_featurecounts_round2.sh
#purpose: mark duplicate reads and quantify STAR-aligned BAM files at the gene level using featureCounts
#To start this job from the sctld_jamboree/tagseq directory, use:
#bsub -P transcriptomics < ./bash/tagseq/unique_featurecounts_2.sh

#BSUB -J feat
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o unique_featurecounts_2_%J.out
#BSUB -e unique_featurecounts_2_%J.err
#BSUB -n 8
#BSUB -u mconnelly@rsmas.miami.edu
#BSUB -N

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/projects/transcriptomics/mikeconnelly"
prodir="/scratch/projects/transcriptomics/mikeconnelly/projects/anti_phys"

# making a list of sample names
samples=$(cat ${prodir}/data/quantseq_samples_2.txt)

module load samtools
module load subread

for samp in $samples ; do
samtools index ${prodir}/outputs/alignments/${samp}_Aligned.sortedByCoord.out.bam
samtools view -q 255 -Sub ${prodir}/outputs/STARalign_Pdam/${samp}_PdamAligned.sortedByCoord.out.bam \
-o ${prodir}/outputs/STARalign_Pdam/${samp}_Aligned.sortedByCoord.out.uniq.bam
done

featureCounts -p -T 8 -t gene \
-g ID \
-a ${mcs}/sequences/genomes/coral/pocillopora/pdam_genome.gff \
-o ${prodir}/outputs/STARcounts_Pdam/PocAnti_Pdam_Round2.counts \
${prodir}/outputs/STARalign_Pdam/*uniq.bam
