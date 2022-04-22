#!/bin/bash
#./bash/tagseq/staralignwrapper_symD.sh
#purpose: align quantseq files against the Cladocopium and Durusdinium concatenated reference genome
#To start this job from the anti_phys directory, use:
#bsub -P transcriptomics < ./bash/tagseq/staralignwrapper_symbio.sh

#BSUB -J staralignwrapper_symbio
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o staralignwrapper_symbio%J.out
#BSUB -e staralignwrapper_symbio%J.err
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
echo '#!/bin/bash' > "${prodir}"/bash/jobs/"${sample}"_staralign_symbio.job
echo '#BSUB -q bigmem' >> "${prodir}"/bash/jobs/"${sample}"_staralign_symbio.job
echo '#BSUB -J '"${sample}"_staralign_symbio'' >> "${prodir}"/bash/jobs/"${sample}"_staralign_symbio.job
echo '#BSUB -o '"${prodir}"/outputs/logfiles/"$sample"staralign_symbio%J.out'' >> "${prodir}"/bash/jobs/"${sample}"_staralign_symbio.job
echo '#BSUB -e '"${prodir}"/outputs/errorfiles/"$sample"staralign_symbio%J.err'' >> "${prodir}"/bash/jobs/"${sample}"_staralign_symbio.job
echo '#BSUB -n 8' >> "${prodir}"/bash/jobs/"${sample}"_staralign_symbio.job
echo '#BSUB -W 4:00' >> "${prodir}"/bash/jobs/"${sample}"_staralign_symbio.job
echo '#BSUB -u mconnelly@rsmas.miami.edu' >> "${prodir}"/bash/jobs/"${sample}"_staralign_symbio.job
echo '#BSUB -N' >> "${prodir}"/bash/jobs/"${sample}"_staralign_symbio.job

#   input command to run STAR aligner
echo ${mcs}/programs/STAR-2.5.3a/bin/Linux_x86_64/STAR \
--runMode alignReads \
--quantMode TranscriptomeSAM \
--runThreadN 16 \
--readFilesIn ${prodir}/outputs/STARalign_Pdam/${sample}_PdamUnmapped.out.mate1 \
--genomeDir ${mcs}/sequences/genomes/symbiodiniaceae/all_sym/STARindex \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile  ${mcs}/sequences/genomes/symbiodiniaceae/all_sym/sym_cat_genome.gff \
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
--outFileNamePrefix ${prodir}/outputs/STARalign_Symbio/${sample}_Symbio >> "${prodir}"/bash/jobs/"${sample}"_staralign_symbio.job

#lets me know file is done
echo 'echo' "STAR alignment of $sample complete" >> "${prodir}"/bash/jobs/"${sample}"_staralign_symbio.job
echo "STAR alignment script of $sample submitted"
#   submit generated trimming script to job queue
bsub < "${prodir}"/bash/jobs/"${sample}"_staralign_symbio.job
done
