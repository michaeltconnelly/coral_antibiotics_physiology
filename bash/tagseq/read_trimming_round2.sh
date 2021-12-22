#!/bin/bash
#./bash/tagseq/read_trimming_round2.sh
#purpose: create wrapper scripts for read trimming using bbduk from BBtools suite
#To start this job from the anti_phys directory, use:
#bsub -P transcriptomics < ./bash/tagseq/read_trimming_round2.sh

#BSUB -J read_trimming_round2
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o read_trimming_round2_%J.out
#BSUB -e read_trimming_round2_%J.err
#BSUB -n 8
#BSUB -u mconnelly@rsmas.miami.edu
#BSUB -N

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/projects/transcriptomics/mikeconnelly"
prodir="/scratch/projects/transcriptomics/mikeconnelly/projects/anti_phys"

# making a list of sample names
samples=$(cat ${prodir}/data/quantseq_samples_2.txt)

#lets me know which files are being processed
echo "These are the samples to be trimmed:"
echo $samples

#loop to automate generation of scripts to direct sequence file trimming
for sample in $samples
do \
echo "Preparing script for ${sample}"

#   input BSUB commands
echo '#!/bin/bash' > ${prodir}/bash/jobs/"${sample}"_trim.job
echo '#BSUB -q general' >> ${prodir}/bash/jobs/"${sample}"_trim.job
echo '#BSUB -J '"${sample}"'_trim' >> ${prodir}/bash/jobs/"${sample}"_trim.job
echo '#BSUB -o "'${prodir}'"/outputs/logfiles/'"${sample}"'_trim.out' >> ${prodir}/bash/jobs/"${sample}"_trim.job
echo '#BSUB -e "'${prodir}'"/outputs/errorfiles/'"${sample}"'_trim.err' >> ${prodir}/bash/jobs/"${sample}"_trim.job

#   input command to load modules for trimming
echo 'module load java/1.8.0_60' >> ${prodir}/bash/jobs/"${sample}"_trim.job

#   input command to unzip raw reads before trimming
echo 'echo 'Unzipping "${sample}"'' >> "${prodir}"/bash/jobs/"${sample}"_trim.job
echo 'gunzip '"${prodir}"/data/quantseq_reads/"${sample}".fastq.gz >> "${prodir}"/bash/jobs/"${sample}"_trim.job

#   input command to trim raw reads
echo 'echo 'Trimming "${sample}"'' >> "${prodir}"/bash/jobs/"${sample}"_trim.job

echo '~/programs/bbmap/bbduk.sh -Xmx512m \
in='"${prodir}"'/data/quantseq_reads/'"${sample}"'_R1.fastq \
out='"${prodir}"'/outputs/trimmed_reads/'"${sample}"'_trimmed.fastq \
ref='"${prodir}"'/data/bbmap_resources/polyA.fa.gz,'"${prodir}"'/data/bbmap_resources/truseq_rna.fa.gz \
k=13 \
ktrim=r \
useshortkmers=T \
mink=5 \
qtrim=10 \
minlength=20' >> ${prodir}/bash/jobs/"${sample}"_trim.job

echo 'echo '"${sample}"' successfully trimmed' >> "${prodir}"/bash/jobs/"${sample}"_trim.job

bsub < ${prodir}/bash/jobs/"${sample}"_trim.job
done
