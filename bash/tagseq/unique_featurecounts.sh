# /bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -pe mthread 16
#$ -q sThC.q
#$ -l mres=64G,h_data=4G,h_vmem=4G
#$ -j y
#$ -N unique_featurecounts
#$ -o /scratch/nmnh_corals/connellym/projects/anti_phys/bash/jobs/unique_featurecounts.log
#$ -m bea
#$ -M connellym@si.edu
#
# ----------------Modules------------------------- #
module load bioinformatics/samtools
#
# ----------------Your Commands------------------- #
#!/bin/bash
#./bash/tagseq/unique_featurecounts.sh
#purpose: mark duplicate reads and quantify STAR-aligned BAM files at the gene level using featureCounts
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
echo + NSLOTS = $NSLOTS

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/anti_phys"

# making a list of sample names
files=$(ls /scratch/nmnh_corals/connellym/projects/anti_phys/data/quantseq_reads/)
samples=$(echo "$files" | cut -d . -f 1 | sort -u)

for sample in $samples ; do
samtools index ${prodir}/outputs/STARalign_Pdam/${sample}_Pdam_Aligned.sortedByCoord.out.bam
samtools view -q 255 -Sub ${prodir}/outputs/STARalign_Pdam/${sample}_Pdam_Aligned.sortedByCoord.out.bam \
-o ${prodir}/outputs/STARalign_Pdam/${sample}_Aligned.sortedByCoord.out.uniq.bam
done

${mcs}/programs/subread-2.0.3-source/bin/featureCounts -p -T 8 -t gene \
-g ID \
-a /home/connellym/sequences/pdam_annotation.gff3 \
-o ${prodir}/outputs/STARcounts_Pdam/PocAnti_Pdam.counts \
${prodir}/outputs/STARalign_Pdam/*uniq.bam
#
echo = `date` job $JOB_NAME done
