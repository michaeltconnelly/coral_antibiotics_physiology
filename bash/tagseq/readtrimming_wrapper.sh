#!/bin/bash
#./bash/tagseq/readtrimming_wrapper.sh
#purpose: create wrapper scripts for read trimming using bbduk from BBtools suite

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/anti_phys"

# making a list of sample names
files=$(ls /scratch/nmnh_corals/connellym/projects/anti_phys/data/quantseq_reads/)
samples=$(echo "$files" | cut -d . -f 1 | sort -u)
#samples=$(cat ${prodir}/data/quantseq_samples.txt ${prodir}/data/quantseq_samples_2.txt)

#lets me know which files are being processed
echo "These are the samples to be trimmed:"
echo $samples

#loop to automate generation of scripts to direct sequence file trimming
for sample in $samples
do \
echo "Preparing script for ${sample}"
#   input QSUB commands
echo "# /bin/sh" > ${prodir}/bash/jobs/${sample}_trim.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/${sample}_trim.job
echo "#$  -S /bin/sh
#$ -pe mthread 16
#$ -q sThC.q
#$ -l mres=64G,h_data=4G,h_vmem=4G" >> ${prodir}/bash/jobs/${sample}_trim.job
echo "#$ -j y
#$ -N trimming_${sample}
#$ -o ${prodir}/bash/jobs/${sample}_trim.log
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/${sample}_trim.job
#
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/${sample}_trim.job
echo "module load bioinformatics/trimmomatic" >> ${prodir}/bash/jobs/${sample}_trim.job
#
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/${sample}_trim.job
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/${sample}_trim.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/${sample}_trim.job
#
#   input command to load modules for trimming
echo 'module load java/1.8' >> ${prodir}/bash/jobs/${sample}_trim.job

#   input command to unzip raw reads before trimming
echo 'echo 'Unzipping ${sample}'' >> ${prodir}/bash/jobs/${sample}_trim.job
echo 'gunzip '${prodir}/data/quantseq_reads/${sample}.fastq.gz >> ${prodir}/bash/jobs/${sample}_trim.job

#   input command to trim raw reads
echo 'echo 'Trimming ${sample}'' >> "${prodir}"/bash/jobs/${sample}_trim.job
#
echo ${mcs}'/programs/bbmap/bbduk.sh -Xmx512m \
in='${prodir}'/data/quantseq_reads/'${sample}'.fastq \
out='${prodir}'/outputs/trimmed_reads/'${sample}'_trimmed.fastq \
ref='${prodir}'/data/bbmap_resources/polyA.fa.gz,'${prodir}'/data/bbmap_resources/truseq_rna.fa.gz \
k=13 \
ktrim=r \
useshortkmers=T \
mink=5 \
qtrim=10 \
minlength=20' >> ${prodir}/bash/jobs/${sample}_trim.job
#
echo 'echo '${sample}' successfully trimmed' >> "${prodir}"/bash/jobs/${sample}_trim.job
#
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/${sample}_trim.job
# submit job
qsub ${prodir}/bash/jobs/${sample}_trim.job
#
done
