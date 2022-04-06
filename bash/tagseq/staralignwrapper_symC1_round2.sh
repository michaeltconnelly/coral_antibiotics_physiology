#!/bin/bash
#./bash/tagseq/staralignwrapper_symC1_round2.sh
#purpose: align quantseq files against the Pocillopora reference genome
#To start this job from the anti_phys directory, use:
#bsub -P transcriptomics < ./bash/tagseq/staralignwrapper_symC1_round2.sh

#BSUB -J staralignwrapper_symC1_round2
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o staralignwrapper_symC1%J.out
#BSUB -e staralignwrapper_symC1%J.err
#BSUB -n 8
#BSUB -u mconnelly@rsmas.miami.edu
#BSUB -N

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/projects/transcriptomics/mikeconnelly"
prodir="/scratch/projects/transcriptomics/mikeconnelly/projects/anti_phys"

# making a list of sample names
samples=$(cat ${prodir}/data/quantseq_samples_2.txt)

#lets me know which files are being processed
echo "These are the reads to be aligned to the Cladocopium goreaui reference genome: $samples"

#loop to automate generation of scripts to direct sequence file trimming
for sample in $samples
do \
echo "Aligning ${sample}"

#   input BSUB commands
echo '#!/bin/bash' > "${prodir}"/bash/jobs/"${sample}"_staralign_symC1.job
echo '#BSUB -q bigmem' >> "${prodir}"/bash/jobs/"${sample}"_staralign_symC1.job
echo '#BSUB -J '"${sample}"_staralign_symC1'' >> "${prodir}"/bash/jobs/"${sample}"_staralign_symC1.job
echo '#BSUB -o '"${prodir}"/outputs/logfiles/"$sample"staralign_symC1%J.out'' >> "${prodir}"/bash/jobs/"${sample}"_staralign_symC1.job
echo '#BSUB -e '"${prodir}"/outputs/errorfiles/"$sample"staralign_symC1%J.err'' >> "${prodir}"/bash/jobs/"${sample}"_staralign_symC1.job
echo '#BSUB -n 8' >> "${prodir}"/bash/jobs/"${sample}"_staralign_symC1.job
echo '#BSUB -W 4:00' >> "${prodir}"/bash/jobs/"${sample}"_staralign_symC1.job

#   input command to run STAR aligner
echo ${mcs}/programs/STAR-2.5.3a/bin/Linux_x86_64/STAR \
--runMode alignReads \
--quantMode TranscriptomeSAM \
--runThreadN 16 \
--readFilesIn ${prodir}/outputs/STARalign_Pdam/${sample}_PdamUnmapped.out.mate1 \
--genomeDir ${mcs}/sequences/genomes/symbiodinium/STARindex \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile  ${mcs}/sequences/genomes/symbiodinium/symC1_genome.gff \
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
--outFileNamePrefix ${prodir}/outputs/STARalign_SymC1/${sample}_SymC1 >> "${prodir}"/bash/jobs/"${sample}"_staralign_symC1.job

#lets me know file is done
echo 'echo' "STAR alignment of $sample complete" >> "${prodir}"/bash/jobs/"${EAPSIsample}"_staralign_symC1.job
echo "STAR alignment script of $sample submitted"
#   submit generated trimming script to job queue
bsub < "${prodir}"/bash/jobs/"${sample}"_staralign_symC1.job
done
