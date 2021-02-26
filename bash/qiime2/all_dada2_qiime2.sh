#!/bin/bash
#purpose: Use DADA2 to dereplicate and merge paired-end reads and create feature table of ASVs
#conda activate qiime2-2019.10

prodir="/Users/mikeconnelly/computing/projects/anti_phys"
#exp="AxH"
echo "DADA2 dereplicaton and summarization process started"

### Denoise, dereplicate and merge paired-end reads using DADA2
### Create feature table and representative sequences
 qiime dada2 denoise-paired \
 --i-demultiplexed-seqs ${prodir}/outputs/qiime2/qza/all.qza \
 --p-trunc-len-f 200 \
 --p-trunc-len-r 180 \
 --p-trim-left-f 5 \
 --p-trim-left-r 5 \
 --p-trunc-q 20 \
 --o-table ${prodir}/outputs/qiime2/qza/all_table.qza \
 --o-representative-sequences ${prodir}/outputs/qiime2/qza/all_rep-seqs.qza \
 --o-denoising-stats ${prodir}/outputs/qiime2/qza/all_denoise-stats.qza \

 qiime feature-table summarize \
 --i-table ${prodir}/outputs/qiime2/qza/all_table.qza \
 --o-visualization ${prodir}/outputs/qiime2/qzv/all_table.qzv \
 --m-sample-metadata-file ${prodir}/data/qiime2_metadata.tsv

 qiime feature-table tabulate-seqs \
 --i-data ${prodir}/outputs/qiime2/qza/all_rep-seqs.qza \
 --o-visualization ${prodir}/outputs/qiime2/qzv/all_rep-seqs.qzv
###      854 representative sequences (features) after denoising and dereplication
 qiime metadata tabulate \
 --m-input-file ${prodir}/outputs/qiime2/qza/all_denoise-stats.qza \
 --o-visualization ${prodir}/outputs/qiime2/qzv/all_denoise-stats.qzv
###      Between 704 - 21,152 non-chimeric reads per sample

bash ${prodir}/bash/qiime2/all_classify_qiime2.sh
