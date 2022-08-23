#!/bin/bash
#./bash/quantseq_setup_start.sh

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/nmnh_corals/connellym"
prodir="/scratch/nmnh_corals/connellym/projects/anti_phys"

samples=$(cat quantseq_samples.txt)
echo "Pipeline setup process started"

#copy program binaries and change permissions
cp -r ~/programs/bbmap ${mcs}/programs
#cp -r ~/programs/FastQC ${mcs}/programs
cp -r ~/programs/STAR-2.7.9a ${mcs}/programs
cp -r ~/programs/subread-2.0.3-source ${mcs}/programs
echo "Program files copied to scratch"

#make file structure for pipeline file input/output
mkdir ${prodir}/bash/jobs
mkdir ${prodir}/data/bbmap_resources
mkdir ${prodir}/outputs/logfiles
mkdir ${prodir}/outputs/errorfiles
mkdir ${prodir}/outputs/fastqcs
mkdir ${prodir}/outputs/trimqcs
mkdir ${prodir}/outputs/trimmed_reads
mkdir ${prodir}/outputs/STARalign_Pdam
mkdir ${prodir}/outputs/STARcounts_Pdam
echo "Filesystem and project directories created"

cp ~/programs/bbmap/resources/polyA.fa.gz ${prodir}/data/bbmap_resources
cp ~/programs/bbmap/resources/truseq_rna.fa.gz ${prodir}/data/bbmap_resources

#Call first scripts in analysis pipeline
#bsub -P transcriptomics < ${prodir}/bash/fastqc.sh
#bsub -P transcriptomics < ${prodir}/bash/trimmomatic.sh
echo "RNAseq pipeline scripts successfully activated"
