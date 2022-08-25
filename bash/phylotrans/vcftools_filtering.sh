# /bin/sh
# ----------------Parameters---------------------- #
#$  -S /bin/sh
#$ -q sThC.q
#$ -l mres=4G,h_data=4G,h_vmem=4G
#$ -cwd
#$ -j y
#$ -N vcftools_filter_primary
#$ -o /scratch/nmnh_corals/connellym/projects/anti_phys/bash/jobs/vcftools_filter_primary.log
#$ -m bea
#$ -M connellym@si.edu
#
# ----------------Modules------------------------- #
module load bioinformatics/vcflib
module load bioinformatics/vcftools
#
# ----------------Your Commands------------------- #
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
VCF_NAME="$1"
prodir="/scratch/nmnh_corals/connellym/projects/anti_phys"
VCF_IN="${prodir}/outputs/phylotrans_Pdam/${VCF_NAME}.vcf.gz"
VCF_OUT="${prodir}/outputs/phylotrans_Pdam/${VCF_NAME}_filtered_primary.vcf.gz"
# set basic filtering parameters
MAF="0.05"
MISS="0.5"
QUAL="30"
MIN_DEPTH="20"
MAX_DEPTH="5000"
THIN_BP="5000"
# run vcftools command
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --thin $THIN_BP --recode --stdout | gzip -c > \
$VCF_OUT
#
echo = `date` job $JOB_NAME done
