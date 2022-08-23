# /bin/sh
# ----------------Parameters---------------------- #
#$  -S /bin/sh
#$ -pe mthread 8
#$ -q mThC.q
#$ -l mres=64G,h_data=8G,h_vmem=8G
#$ -cwd
#$ -j y
#$ -N gatk_genotypegvcfs_all
#$ -o gatk_genotypegvcfs_all.log
#$ -m bea
#$ -M connellym@si.edu
#
# ----------------Modules------------------------- #
module load bioinformatics/gatk
#
# ----------------Your Commands------------------- #
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
echo + NSLOTS = $NSLOTS
#
prodir="/scratch/nmnh_corals/connellym/projects/anti_phys"
dir="/scratch/nmnh_corals/connellym/projects/anti_phys/outputs/phylotrans_Pdam"
#
java -jar /share/apps/bioinformatics/gatk/3.8.1.0/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R /home/connellym/sequences/pdam_scaffolds.fasta \
-V ${prodir}/outputs/phylotrans_Pdam/samples_all.g.vcf.gz \
-o ${prodir}/outputs/phylotrans_Pdam/samples_all.vcf.gz
#
echo = `date` job $JOB_NAME done
