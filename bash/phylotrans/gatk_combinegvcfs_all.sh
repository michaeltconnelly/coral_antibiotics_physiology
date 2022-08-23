# /bin/sh
# ----------------Parameters---------------------- #
#$  -S /bin/sh
#$ -pe mthread 16
#$ -q mThC.q
#$ -l mres=64G,h_data=8G,h_vmem=8G
#$ -cwd
#$ -j y
#$ -N gatk_combinegvcfs_all
#$ -o /scratch/nmnh_corals/connellym/projects/anti_phys/bash/jobs/gatk_combinegvcfs_all.log
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
--variant ${dir}/PAN-10_Baseline_007_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-10_Control_061_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-10_Control_062_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-10_Control_064_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-32_Antibiotics_121_Pdam.g.vcf.gz \
--variant ${dir}/PAN-32_Antibiotics_122_Pdam.g.vcf.gz \
--variant ${dir}/PAN-32_Antibiotics_124_Pdam.g.vcf.gz \
--variant ${dir}/PAN-32_Baseline_009_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-32_Baseline_011_Pdam.g.vcf.gz \
--variant ${dir}/PAN-32_Baseline_012_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-32_Control_066_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-32_Control_067_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-32_Control_068_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-34_Antibiotics_127_Pdam.g.vcf.gz \
--variant ${dir}/PAN-34_Baseline_016_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-34_Control_071_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-35_Antibiotics_129_Pdam.g.vcf.gz \
--variant ${dir}/PAN-35_Antibiotics_130_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-35_Antibiotics_131_Pdam.g.vcf.gz \
--variant ${dir}/PAN-35_Antibiotics_132_Pdam.g.vcf.gz \
--variant ${dir}/PAN-35_Baseline_017_Pdam.g.vcf.gz \
--variant ${dir}/PAN-35_Baseline_018_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-35_Baseline_020_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-35_Control_073_Pdam.g.vcf.gz \
--variant ${dir}/PAN-35_Control_074_Pdam.g.vcf.gz \
--variant ${dir}/PAN-35_Control_075_Pdam.g.vcf.gz \
--variant ${dir}/PAN-35_Control_076_Pdam.g.vcf.gz \
--variant ${dir}/PAN-37_Antibiotics_133_Pdam.g.vcf.gz \
--variant ${dir}/PAN-37_Antibiotics_134_Pdam.g.vcf.gz \
--variant ${dir}/PAN-37_Antibiotics_135_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-37_Antibiotics_136_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-37_Baseline_021_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-37_Baseline_023_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-37_Control_078_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-39_Antibiotics_137_Pdam.g.vcf.gz \
--variant ${dir}/PAN-39_Antibiotics_138_Pdam.g.vcf.gz \
--variant ${dir}/PAN-39_Antibiotics_139_Pdam.g.vcf.gz \
--variant ${dir}/PAN-39_Antibiotics_140_Pdam.g.vcf.gz \
--variant ${dir}/PAN-39_Baseline_025_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-39_Baseline_026_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-39_Baseline_027_Pdam.g.vcf.gz \
--variant ${dir}/PAN-39_Baseline_028_Pdam.g.vcf.gz \
--variant ${dir}/PAN-39_Control_081_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-39_Control_082_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-39_Control_083_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-39_Control_084_Pdam.g.vcf.gz \
--variant ${dir}/PAN-41_Antibiotics_141_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-41_Antibiotics_144_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-41_Baseline_031_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-41_Baseline_032_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-43_Antibiotics_145_Pdam.g.vcf.gz \
--variant ${dir}/PAN-43_Antibiotics_146_Pdam.g.vcf.gz \
--variant ${dir}/PAN-43_Antibiotics_147_Pdam.g.vcf.gz \
--variant ${dir}/PAN-43_Antibiotics_148_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-43_Baseline_034_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-44_Antibiotics_149_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-44_Antibiotics_151_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-44_Antibiotics_152_Pdam.g.vcf.gz \
--variant ${dir}/PAN-44_Baseline_037_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-78_Antibiotics_157_Pdam.g.vcf.gz \
--variant ${dir}/PAN-78_Antibiotics_158_Pdam.g.vcf.gz \
--variant ${dir}/PAN-78_Antibiotics_159_Pdam.g.vcf.gz \
--variant ${dir}/PAN-78_Antibiotics_160_Pdam.g.vcf.gz \
--variant ${dir}/PAN-78_Baseline_045_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-78_Baseline_046_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-78_Baseline_048_Pdam.g.vcf.gz \
--variant ${dir}/PAN-78_Control_101_Pdam.g.vcf.gz \
--variant ${dir}/PAN-78_Control_102_Pdam.g.vcf.gz \
--variant ${dir}/PAN-78_Control_103_Pdam.g.vcf.gz \
--variant ${dir}/PAN-78_Control_104_Pdam.g.vcf.gz \
--variant ${dir}/PAN-79_Antibiotics_161_Pdam.g.vcf.gz \
--variant ${dir}/PAN-79_Antibiotics_162_Pdam.g.vcf.gz \
--variant ${dir}/PAN-79_Antibiotics_163_Pdam.g.vcf.gz \
--variant ${dir}/PAN-79_Antibiotics_164_Pdam.g.vcf.gz \
--variant ${dir}/PAN-79_Baseline_049_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-79_Baseline_050_Pdam.g.vcf.gz \
--variant ${dir}/PAN-79_Baseline_051_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-79_Baseline_052_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-79_Control_105_Pdam.g.vcf.gz \
--variant ${dir}/PAN-79_Control_106_Pdam.g.vcf.gz \
--variant ${dir}/PAN-79_Control_107_Pdam.g.vcf.gz \
--variant ${dir}/PAN-79_Control_108_Pdam.g.vcf.gz \
--variant ${dir}/PAN-83_Antibiotics_165_Pdam.g.vcf.gz \
--variant ${dir}/PAN-83_Antibiotics_166_Pdam.g.vcf.gz \
--variant ${dir}/PAN-83_Antibiotics_167_Pdam.g.vcf.gz \
--variant ${dir}/PAN-83_Antibiotics_168_Pdam.g.vcf.gz \
--variant ${dir}/PAN-83_Baseline_054_R1_Pdam.g.vcf.gz \
--variant ${dir}/PAN-83_Baseline_056_Pdam.g.vcf.gz \
--variant ${dir}/PAN-83_Control_109_Pdam.g.vcf.gz \
--variant ${dir}/PAN-83_Control_110_Pdam.g.vcf.gz \
--variant ${dir}/PAN-83_Control_111_Pdam.g.vcf.gz \
--variant ${dir}/PAN-83_Control_112_Pdam.g.vcf.gz \
--variant ${dir}/URA-51_Antibiotics_153_Pdam.g.vcf.gz \
--variant ${dir}/URA-51_Antibiotics_154_R1_Pdam.g.vcf.gz \
--variant ${dir}/URA-51_Antibiotics_155_Pdam.g.vcf.gz \
--variant ${dir}/URA-51_Antibiotics_156_Pdam.g.vcf.gz \
--variant ${dir}/URA-51_Control_097_Pdam.g.vcf.gz \
--variant ${dir}/URA-51_Control_098_Pdam.g.vcf.gz \
--variant ${dir}/URA-51_Control_099_Pdam.g.vcf.gz \
--variant ${dir}/URA-51_Control_100_Pdam.g.vcf.gz \ \
--out ${dir}/outputs/phylotrans_Pdam/samples_all.g.vcf.gz
#
echo = `date` job $JOB_NAME done
#
qsub ${prodir}/jobs/gatk_genotypegvcfs_all.sh
