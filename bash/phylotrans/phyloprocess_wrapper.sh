#!/bin/bash
#./bash/tagseq/phyloprocess_wrapper.sh
#purpose: index bam files, mark duplicate reads and splitN for downstream SNP calling

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/anti_phys"

# making a list of sample names
files=$(ls /scratch/nmnh_corals/connellym/projects/anti_phys/data/quantseq_reads/)
samples=$(echo "$files" | cut -d . -f 1 | sort -u)

for sample in $samples
do
echo "# /bin/sh" > ${prodir}/bash/jobs/phyloprocess_${sample}.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
echo "#$  -S /bin/sh
#$ -pe mthread 16
#$ -q sThC.q
#$ -l mres=32G,h_data=2G,h_vmem=2G
#$ -cwd
#$ -j y
#$ -N phyloprocess_${sample}
#$ -o ${prodir}/bash/jobs/${sample}_phyloprocess.log
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
#
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
echo "module load bioinformatics/picard-tools" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
echo "module load bioinformatics/samtools" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
echo "module load bioinformatics/gatk" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
#
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
echo "#" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
#
echo 'echo "This is the sample being processed:"' >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
echo "echo $sample" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
echo 'echo "Starting samtools sort and index steps"' >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
#
echo "samtools sort \
${prodir}/outputs/STARalign_Pdam/${sample}_Pdam_Aligned.out.bam \
> ${prodir}/outputs/phylotrans_Pdam/${sample}_Pdam_Aligned.sorted.out.bam" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
echo "#" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
#
echo "samtools index -b \
${prodir}/outputs/phylotrans_Pdam/${sample}_Pdam_Aligned.sorted.out.bam" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
echo "#" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
#
echo 'echo "Starting PicardTools and GATK processing steps"' >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
echo "java -jar /share/apps/bioinformatics/picard-tools/2.20.6/picard.jar \
AddOrReplaceReadGroups \
INPUT=${prodir}/outputs/phylotrans_Pdam/${sample}_Pdam_Aligned.sorted.out.bam \
OUTPUT=${prodir}/outputs/phylotrans_Pdam/${sample}_Pdam_Aligned.sorted.out.rg.bam \
RGID=id \
RGLB=library \
RGPL=illumina \
RGPU=unit1 \
RGSM=${sample}" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
echo "#" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
#
echo "java -jar /share/apps/bioinformatics/picard-tools/2.20.6/picard.jar \
MarkDuplicates \
INPUT=${prodir}/outputs/phylotrans_Pdam/${sample}_Pdam_Aligned.sorted.out.rg.bam \
OUTPUT=${prodir}/outputs/phylotrans_Pdam/${sample}_Pdam_Aligned.sorted.out.md.rg.bam \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
METRICS_FILE=${prodir}/outputs/phylotrans_Pdam/${sample}_marked_dup_metrics.txt" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
echo "#" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
#
echo "java -jar /share/apps/bioinformatics/gatk/3.8.1.0/GenomeAnalysisTK.jar \
-T SplitNCigarReads \
-I ${prodir}/outputs/phylotrans_Pdam/${sample}_Pdam_Aligned.sorted.out.md.rg.bam \
-o ${prodir}/outputs/phylotrans_Pdam/${sample}_Pdam_Aligned.sorted.out.md.rg.splitN.bam \
-R /home/connellym/sequences/pdam_scaffolds.fasta \
-rf ReassignOneMappingQuality \
-RMQF 255 -RMQT 60 \
-U ALLOW_N_CIGAR_READS" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
echo "#" >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/phyloprocess_${sample}.job
qsub ${prodir}/bash/jobs/phyloprocess_${sample}.job
done
