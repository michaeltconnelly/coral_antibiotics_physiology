#!/bin/bash
#./bash/tagseq/unique_featurecounts_symbio.sh
#purpose: mark duplicate reads and quantify STAR-aligned BAM files at the gene level using featureCounts
#To start this job from the anti_phys directory, use:
#bsub -P transcriptomics < ./bash/tagseq/unique_featurecounts_symbio.sh

#BSUB -J unique_featurecounts_symbio
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o unique_featurecounts_symbio%J.out
#BSUB -e unique_featurecounts_symbio%J.err
#BSUB -n 8
#BSUB -u mconnelly@rsmas.miami.edu
#BSUB -N

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/projects/transcriptomics/mikeconnelly"
prodir="/scratch/projects/transcriptomics/mikeconnelly/projects/anti_phys"

# making a list of sample names
samples=$(cat ${prodir}/data/quantseq_samples.txt ${prodir}/data/quantseq_samples_2.txt)

module load samtools
module load subread

for sample in $samples ; do
samtools index ${prodir}/outputs/STARalign_Symbio/${sample}_SymbioAligned.sortedByCoord.out.bam
samtools view -q 255 -Sub ${prodir}/outputs/STARalign_Symbio/${sample}_SymbioAligned.sortedByCoord.out.bam \
-o ${prodir}/outputs/STARalign_Symbio/${sample}_SymbioAligned.sortedByCoord.out.uniq.bam
echo "${sample} processed with samtools index, view"
done

featureCounts -p -T 8 -t gene \
-g ID \
-a ${mcs}/sequences/genomes/symbiodiniaceae/all_sym/sym_cat_genome.gff \
-o ${prodir}/outputs/STARcounts_Symbio/PocAnti_Symbio.counts \
${prodir}/outputs/STARalign_Symbio/*uniq.bam
