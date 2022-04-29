#!/bin/bash
#./bash/tagseq/symalignwrapper.sh
#purpose: align quantseq files against the Symbiodinium, Breviolum, Cladocopium and Durusdinium concatenated reference genome
#To start this job from the anti_phys directory, use:
#bsub -P transcriptomics < ./bash/tagseq/staralignwrapper_symbio.sh

#BSUB -J symalignwrapper
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o symalignwrapper_%J.out
#BSUB -e symalignwrapper_%J.err
#BSUB -n 8
#BSUB -u mconnelly@rsmas.miami.edu
#BSUB -N

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/projects/transcriptomics/mikeconnelly"
prodir="/scratch/projects/transcriptomics/mikeconnelly/projects/anti_phys"

# making a list of sample names
samples=$(cat ${prodir}/data/quantseq_samples.txt ${prodir}/data/quantseq_samples_2.txt)

#lets me know which files are being processed
echo "These are the reads to be aligned to the concatenated Symbiodiniaceae reference genome: $samples"

#loop to automate generation of scripts to direct sequence file trimming
for sample in $samples
do \
echo "Aligning ${sample}"

#   input BSUB commands
echo '#!/bin/bash' > "${prodir}"/bash/jobs/"${sample}"_symalign_ABCD.job
echo '#BSUB -q general' >> "${prodir}"/bash/jobs/"${sample}"_symalign_ABCD.job
echo '#BSUB -J '"${sample}"_staralign_symbio'' >> "${prodir}"/bash/jobs/"${sample}"_symalign_ABCD.job
echo '#BSUB -o '"${prodir}"/outputs/logfiles/"$sample"staralign_symbio%J.out'' >> "${prodir}"/bash/jobs/"${sample}"_symalign_ABCD.job
echo '#BSUB -e '"${prodir}"/outputs/errorfiles/"$sample"staralign_symbio%J.err'' >> "${prodir}"/bash/jobs/"${sample}"_symalign_ABCD.job
echo '#BSUB -n 8' >> "${prodir}"/bash/jobs/"${sample}"_symalign_ABCD.job
echo '#BSUB -W 4:00' >> "${prodir}"/bash/jobs/"${sample}"_symalign_ABCD.job
echo '#BSUB -u mconnelly@rsmas.miami.edu' >> "${prodir}"/bash/jobs/"${sample}"_symalign_ABCD.job
echo '#BSUB -N' >> "${prodir}"/bash/jobs/"${sample}"_symalign_ABCD.job

echo "module load bowtie2/2.2.6"  >> "${prodir}"/bash/jobs/"${sample}"_symalign_ABCD.job
#   input command to run STAR aligner
echo "bowtie2 \
-x ${mcs}/sequences/genomes/symbiodiniaceae/all_sym/symbiosymABCDgenome \
-U ${prodir}/outputs/STARalign_Pdam/${sample}_PdamUnmapped.out.mate1 \
--local \
-p 8 \
-S ${prodir}/outputs/symABCD/${sample}_symABCD.sam" >> "${prodir}"/bash/jobs/"${sample}"_symalign_ABCD.job

#lets me know file is done
echo 'echo' "Bowtie alignment of $sample complete" >> "${prodir}"/bash/jobs/"${sample}"_symalign_ABCD.job
echo "Bowtie alignment script of $sample submitted"
#   submit generated trimming script to job queue
bsub < "${prodir}"/bash/jobs/"${sample}"_symalign_ABCD.job
done
