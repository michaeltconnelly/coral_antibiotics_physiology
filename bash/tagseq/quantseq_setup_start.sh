#!/bin/bash
#./bash/quantseq_setup_start.sh
#purpose: create file directory structure in Pegasus scratch space and copy key reference files, start pipeline
#To start this job from the anti_phys directory, use:
#bsub -P transcriptomics < ./bash/quantseq_setup_start.sh

#BSUB -J quantseq_setup_start
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o quantseq_setup_start%J.out
#BSUB -e quantseq_setup_start%J.err
#BSUB -n 8
#BSUB -u mconnelly@rsmas.miami.edu
#BSUB -N

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/projects/transcriptomics/mikeconnelly"
prodir="/scratch/projects/transcriptomics/mikeconnelly/projects/anti_phys"

samples=$(cat quantseq_samples.txt)
echo "Pipeline setup process started"

#copy reference sequences
#Pocillopora genome - .fasta, .gff, STAR index
cp -r ~/sequences/genomes/coral/pocillopora ${mcs}/sequences/genomes/coral/
#Symbiodinium C1 genome - .fasta, .gff, STAR index
cp -r ~/sequences/genomes/symbiodinium ${mcs}/sequences/genomes/
echo "Reference genome sequences copied to scratch"

#copy program binaries, change permissions, and load necessary modules
#execute FASTQC and Trimmomatic using Pegasus modules
module load java/1.8.0_60
module load trimmomatic/0.36
cp -r ~/programs/bbmap ${mcs}/programs
#cp -r ~/programs/FastQC ${mcs}/programs
cp -r ~/programs/STAR-2.5.3a ${mcs}/programs
cp -r ~/programs/subread-1.6.0-Linux-x86_64 ${mcs}/programs
chmod 755 ${mcs}/programs/FastQC/fastqc
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
bsub -P transcriptomics < ${prodir}/bash/fastqc.sh
bsub -P transcriptomics < ${prodir}/bash/trimmomatic.sh
echo "RNAseq pipeline scripts successfully activated"
