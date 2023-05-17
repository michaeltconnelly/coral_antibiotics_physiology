#!/bin/bash
#./bash/tagseq/phyloprocess_wrapper.sh
#purpose: convert gVCF files to VCF files for use with ANGSD

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/anti_phys"

# making a list of sample names
files=$(ls /scratch/nmnh_corals/connellym/projects/anti_phys/data/quantseq_reads/)
samples=$(echo "$files" | cut -d . -f 1 | sort -u)
# | cut -d _ -f 1,2,3 

for sample in $samples
do
echo "# /bin/sh" > ${prodir}/bash/jobs/gvcf2vcf_${sample}.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/gvcf2vcf_${sample}.job
echo "#$  -S /bin/sh
#$ -pe mthread 16
#$ -q sThC.q
#$ -l mres=32G,h_data=2G,h_vmem=2G
#$ -cwd
#$ -j y
#$ -N gvcf2vcf_${sample}
#$ -o ${prodir}/bash/jobs/${sample}_gvcf2vcf.log
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/gvcf2vcf_${sample}.job
#
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/gvcf2vcf_${sample}.job
echo "module load bioinformatics/bcftools" >> ${prodir}/bash/jobs/gvcf2vcf_${sample}.job
#
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/gvcf2vcf_${sample}.job
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/gvcf2vcf_${sample}.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/gvcf2vcf_${sample}.job
echo "#" >> ${prodir}/bash/jobs/gvcf2vcf_${sample}.job
#
echo 'echo "This is the sample being processed:"' >> ${prodir}/bash/jobs/gvcf2vcf_${sample}.job
echo "echo $sample" >> ${prodir}/bash/jobs/gvcf2vcf_${sample}.job
echo 'echo "Starting conversion"' >> ${prodir}/bash/jobs/gvcf2vcf_${sample}.job
#
echo "bcftools convert --gvcf2vcf \
-f /scratch/nmnh_corals/connellym/sequences/pdam/pdam_genome.fasta \
${prodir}/outputs/phylotrans_Pdam/${sample}_Pdam.g.vcf.gz \
-o ${prodir}/outputs/phylotrans_Pdam/${sample}_Pdam.vcf.gz \
--threads 16" >> ${prodir}/bash/jobs/gvcf2vcf_${sample}.job
echo "#" >> ${prodir}/bash/jobs/gvcf2vcf_${sample}.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/gvcf2vcf_${sample}.job
qsub ${prodir}/bash/jobs/gvcf2vcf_${sample}.job
done