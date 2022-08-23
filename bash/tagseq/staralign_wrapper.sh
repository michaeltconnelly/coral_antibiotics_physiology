#!/bin/bash
#./bash/tagseq/staralign_wrapper.sh
#purpose: align quantseq files against the Pocillopora reference genome

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/anti_phys"

# making a list of sample names
files=$(ls /scratch/nmnh_corals/connellym/projects/anti_phys/data/quantseq_reads/)
samples=$(echo "$files" | cut -d . -f 1 | sort -u)

#lets me know which files are being processed
echo "These are the reads to be aligned to the Pocillopora reference genome:"
echo $samples

#loop to automate generation of scripts to direct sequence file trimming
for sample in $samples
do \
echo "Aligning ${sample}"
echo "# /bin/sh" > ${prodir}/bash/jobs/${sample}_staralign_pdam.job
echo "# ----------------Parameters---------------------- #" >> ${prodir}/bash/jobs/${sample}_staralign_pdam.job
echo "#$  -S /bin/sh
#$ -pe mthread 16
#$ -q sThC.q
#$ -l mres=64G,h_data=4G,h_vmem=4G" >> ${prodir}/bash/jobs/${sample}_staralign_pdam.job
echo "#$ -j y
#$ -N staralign_${sample}
#$ -o ${prodir}/bash/jobs/${sample}_staralign_pdam.log
#$ -m bea
#$ -M connellym@si.edu" >> ${prodir}/bash/jobs/${sample}_staralign_pdam.job
#
echo "# ----------------Modules------------------------- #" >> ${prodir}/bash/jobs/${sample}_staralign_pdam.job
#
echo "# ----------------Your Commands------------------- #" >> ${prodir}/bash/jobs/${sample}_staralign_pdam.job
#
echo 'echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME' >> ${prodir}/bash/jobs/${sample}_staralign_pdam.job
echo 'echo + NSLOTS = $NSLOTS' >> ${prodir}/bash/jobs/${sample}_staralign_pdam.job
#
#   input command to run STAR aligner
echo "${mcs}/programs/STAR-2.7.9a/bin/Linux_x86_64/STAR \
--runMode alignReads \
--quantMode TranscriptomeSAM \
--runThreadN 16 \
--readFilesIn ${prodir}/outputs/trimmed_reads/${sample}_trimmed.fastq \
--genomeDir /home/connellym/sequences/STARidx_pdam \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile /home/connellym/sequences/pdam_annotation.gff3 \
--twopassMode Basic \
--twopass1readsN -1 \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outSAMattributes NH HI NM MD \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix ${prodir}/outputs/STARalign_Pdam/${sample}_Pdam" >> ${prodir}/bash/jobs/${sample}_staralign_pdam.job
#
#lets me know file is done
echo 'echo' "STAR alignment of $sample complete"'' >> ${prodir}/bash/jobs/${sample}_staralign_pdam.job
echo 'echo = `date` job $JOB_NAME done' >> ${prodir}/bash/jobs/${sample}_staralign_pdam.job
#   submit generated trimming script to job queue
qsub ${prodir}/bash/jobs/${sample}_staralign_pdam.job
done
