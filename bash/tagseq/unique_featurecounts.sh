#!/bin/bash
#./bash/tagseq/unique_featurecounts.sh
#purpose: mark duplicate reads and quantify STAR-aligned BAM files at the gene level using featureCounts

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/anti_phys"

# making a list of sample names
files=$(ls /scratch/nmnh_corals/connellym/projects/anti_phys/data/quantseq_reads/)
samples=$(echo "$files" | cut -d . -f 1 | sort -u)

module load bioinformatics/samtools

for sample in $samples ; do
samtools index ${prodir}/outputs/alignments/${sample}_Aligned.sortedByCoord.out.bam
samtools view -q 255 -Sub ${prodir}/outputs/STARalign_Pdam/${samp}_PdamAligned.sortedByCoord.out.bam \
-o ${prodir}/outputs/STARalign_Pdam/${sample}_Aligned.sortedByCoord.out.uniq.bam
done

${mcs}/programs/subread-2.0.3-source/bin/featureCounts -p -T 8 -t gene \
-g ID \
-a /home/connellym/sequences/pdam_annotation.gff3 \
-o ${prodir}/outputs/STARcounts_Pdam/PocAnti_Pdam.counts \
${prodir}/outputs/STARalign_Pdam/*uniq.bam
