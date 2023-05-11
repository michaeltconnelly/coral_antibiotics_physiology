# *Pocillopora* coral antibiotics physiology

These scripts and data are for the analysis of physiology (IPAM fluorometry and respirometry), transcriptome and microbiome data for the publication "**Antibiotics reduce *Pocillopora* coral-associated bacteria diversity, decrease holobiont oxygen consumption, and activate immune gene expression**".

Bacteria 16S and RNAseq data are deposited in the NBCI SRA under BioProject [PRJNA818888](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA818888/). Scripts used for processing bacteria 16S reads with QIIME2 are located in `bash/qiime2` , and scripts used for processing 3'-tagged RNAseq reads are located in `bash/tagseq`. RMarkdown files used in downstream data analysis to produce the manuscript figures are located in `Rmd`.

The following RMarkdown files were used for the manuscript analyses:

-   `mtORF_haplotypes.Rmd`  - Analysis of mtORF sequences to determine *Pocillopora* host lineages

-   `IPAM_analysis.Rmd`  - Analysis of IPAM photosynthetic efficiency data

-   `respR_analysis.Rmd`  - Analysis of microplate respirometry data with respR

-   `phyloseq_microbiome_analysis.Rmd`  - Analysis of bacteria community diversity and composition

-   `quantseq_transcriptome_analysis.Rmd`  - Analysis of *Pocillopora* coral gene expression

![Fig. 1: (a) Map of Panama depicting the collection locations of the fourteen Pocillopora genotypes used in the experiment, and the proportions of different mtORF types from each location. (b) Haplotype networks of the mtORF region for the Pocillopora genotypes used in the antibiotics experiment. Five unique haplotypes (I -- V) were detected in an 828-bp alignment, corresponding to Pocillopora Type 1 (orange) and Type 3 (cyan) sensu Pinzón and LaJeunesse (2011). (c) Representative images of Pocillopora microfragments from each genotype on the first day of the experiment, outlines are colored according to mtORF type, and dominant symbiont type (C/D) is indicated in parentheses. All genotypes in the first row of images were used in transcriptome analyses.](./Fig1.png)
