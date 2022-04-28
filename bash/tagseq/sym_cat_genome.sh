#!/bin/bash
#
# Rename original genome fasta and GFF files after downloading from source

# Create empty files to hold concatenated reference
touch sym_cat_genome.fasta sym_cat_genome.gff
# Append each reference genome into a concatenated reference
cat ../symbiodinium/symA_genome.fasta \
../breviolum/symB_genome.fasta  \
../cladocopium/symC1_genome.fasta \
../durusdinium/symD_genome.fasta >> sym_cat_genome.fasta
# Append each gff files into a concatenated reference
cat ../symbiodinium/symA_genome.gff \
../breviolum/symB_genome.gff \
../cladocopium/symC1_genome.gff \
../durusdinium/symD_genome.gff >> sym_cat_genome.gff
# Remove all but one fasta header for each reference genome, discard all newlines
# chr1, chr2, chr3, chr4
# A,    B,    C,    D,
grep \>.* sym_cat_genome.fasta
#>Smic.scaffold1|size3144590
#scaffold
#SymbC1
#sc0
# find and replace first instance with >chr1, etc.
sed '0,/\>Smic.*/{s/>Smic.*/\>chr1/}' sym_cat_genome.fasta > sym_cat_header1.fasta
#
