#!/bin/bash
#./bash/phylotrans/phylotrans_process_wrapper.sh
#purpose: process aligned Pocillopora reads for transcriptome variant calling
#To start this job from the anti_phys directory, use:
#bsub -P transcriptomics < ./bash/phylotrans/phylotrans_process_wrapper.sh

#BSUB -J phylotrans_process_wrapper
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o phylotrans_process_wrapper_%J.out
#BSUB -e phylotrans_process_wrapper_%J.err
#BSUB -n 8
#BSUB -u mconnelly@rsmas.miami.edu
#BSUB -N

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/projects/transcriptomics/mikeconnelly"
prodir="/scratch/projects/transcriptomics/mikeconnelly/projects/anti_phys"

# making a list of sample names
samples=$(cat ${prodir}/data/quantseq_samples.txt ${prodir}/data/quantseq_samples_2.txt)

# for loop to produce and submit jobs

for sample in $samples
do
#  input BSUB commands
echo "#!/bin/bash
#BSUB -q general
#BSUB -J ${sample}_phylotrans_process_pdam
#BSUB -o ${prodir}/outputs/logfiles/${sample}_phylotrans_process_pdam_%J.out
#BSUB -e ${prodir}/outputs/errorfiles/${sample}_phylotrans_process_%J.err
#BSUB -n 8
#BSUB -W 4:00" >> $prodir/bash/jobs/phylotrans_process_${sample}.job
#
echo "# ----------------Modules------------------------- #" >> $prodir/bash/jobs/phylotrans_process_${sample}.job
echo "module load picard-tools/1.103" >> $prodir/bash/jobs/phylotrans_process_${sample}.job
echo "module load samtools/1.3" >> $prodir/bash/jobs/phylotrans_process_${sample}.job
echo "module load GATK/3.4.0" >> $prodir/bash/jobs/phylotrans_process_${sample}.job

# input command to process files
echo 'echo "This is the sample being processed:"' >> $prodir/bash/jobs/phylotrans_process_${sample}.job
echo "echo $sample" >> $prodir/bash/jobs/phylotrans_process_${sample}.job
echo 'echo "Starting samtools sort and index steps"' >> $prodir/bash/jobs/phylotrans_process_${sample}.job
#
echo "samtools sort \
${prodir}/outputs/STARalign_Pdam/${sample}_PdamAligned.out.bam \
> ${prodir}/outputs/phylotrans_Pdam/${sample}_PdamAligned.sorted.out.bam" >> $prodir/bash/jobs/phylotrans_process_${sample}.job
echo "#" >> $prodir/bash/jobs/phylotrans_process_${sample}.job
#
echo "samtools index -b \
${prodir}/outputs/phylotrans_Pdam/${sample}_PdamAligned.sorted.out.bam" >> $prodir/bash/jobs/phylotrans_process_${sample}.job
echo "#" >> $prodir/bash/jobs/phylotrans_process_${sample}.job
#
echo 'echo "Starting PicardTools and GATK processing steps"' >> $prodir/bash/jobs/phylotrans_process_${sample}.job
echo "java -jar /share/apps/bioinformatics/picard-tools/2.20.6/picard.jar \
AddOrReplaceReadGroups \
INPUT=${prodir}/outputs/phylotrans_Pdam/${sample}_PdamAligned.sorted.out.bam \
OUTPUT=${prodir}/outputs/phylotrans_Pdam/${sample}_PdamAligned.sorted.out.rg.bam \
RGID=id \
RGLB=library \
RGPL=illumina \
RGPU=unit1 \
RGSM=${sample}" >> $prodir/bash/jobs/phylotrans_process_${sample}.job
echo "#" >> $prodir/bash/jobs/phylotrans_process_${sample}.job
#
echo "java -jar /share/apps/bioinformatics/picard-tools/2.20.6/picard.jar \
MarkDuplicates \
INPUT=${prodir}/outputs/phylotrans_Pdam/${sample}_PdamAligned.sorted.out.rg.bam \
OUTPUT=${prodir}/outputs/phylotrans_Pdam/${sample}_PdamAligned.sorted.out.md.rg.bam \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
METRICS_FILE=${prodir}/outputs/phylotrans_Pdam/${sample}_marked_dup_metrics.txt" >> $prodir/bash/jobs/phylotrans_process_${sample}.job
echo "#" >> $prodir/bash/jobs/phylotrans_process_${sample}.job
#
echo "java -jar /share/apps/bioinformatics/gatk/3.8.1.0/GenomeAnalysisTK.jar \
-T SplitNCigarReads \
-I ${prodir}/outputs/phylotrans_Pdam/${sample}_PdamAligned.sorted.out.md.rg.bam \
-o ${prodir}/outputs/phylotrans_Pdam/${sample}_PdamAligned.sorted.out.md.rg.splitN.bam \
-R /scratch/projects/transcriptomics/mikeconnelly/sequences/genomes/coral/pocillopora/pdam_genome.fasta \
-rf ReassignOneMappingQuality \
-RMQF 255 -RMQT 60 \
-U ALLOW_N_CIGAR_READS" >> $prodir/bash/jobs/phylotrans_process_${sample}.job
echo "#" >> $prodir/bash/jobs/phylotrans_process_${sample}.job
#
echo 'echo = `date` job $JOB_NAME done' >> $prodir/bash/jobs/phylotrans_process_${sample}.job
qsub $prodir/bash/jobs/phylotrans_process_${sample}.job
done
