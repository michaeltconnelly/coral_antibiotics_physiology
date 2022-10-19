# /bin/sh
# ----------------Parameters---------------------- #
#$  -S /bin/sh
#$ -pe mthread 16
#$ -q mThC.q
#$ -l mres=64G,h_data=8G,h_vmem=8G
#$ -cwd
#$ -j y
#$ -N gatk_combinegvcfs_genotypes
#$ -o /scratch/nmnh_corals/connellym/projects/anti_phys/bash/jobs/gatk_combinegvcfs_genotypes.log
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
java -jar /share/apps/bioinformatics/gatk/3.8.1.0/GenomeAnalysisTK.jar \
-T CombineGVCFs \
-R /home/connellym/sequences/pdam_scaffolds.fasta \
--variant ${dir}/PAN-05_Baseline_001_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-10_Control_062_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-32_Control_066_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-34_Control_071_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-35_Control_073_Pdam.g.vcf.gz \
--variant ${dir}/PAN-37_Baseline_023_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-39_Control_084_Pdam.g.vcf.gz \
--variant ${dir}/PAN-41_Baseline_031_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-43_Baseline_034_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-44_Baseline_037_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-78_Control_104_Pdam.g.vcf.gz \
--variant ${dir}/PAN-79_Control_108_Pdam.g.vcf.gz \
--variant ${dir}/PAN-83_Control_109_Pdam.g.vcf.gz \
--variant ${dir}/URA-51_Control_097_Pdam.g.vcf.gz \
--out ${prodir}/outputs/phylotrans_Pdam/samples_genotypes.g.vcf.gz
#
echo = `date` job $JOB_NAME done
#
qsub ${prodir}/jobs/gatk_genotypegvcfs_genotypes.sh
