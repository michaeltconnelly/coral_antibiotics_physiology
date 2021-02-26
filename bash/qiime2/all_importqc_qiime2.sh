#!/bin/bash
#purpose: import demultiplexed 16S sequence reads, summarize read depth and quality scores
#conda activate qiime2-2019.10

prodir="/Users/mikeconnelly/computing/projects/anti_phys"
exp="anti_phys"
echo "Import and QC summarize process started"

### Import demultiplexed 16S sequence reads into QIIME2 environment
qiime tools import \
--input-path ${prodir}/data/qiime2_manifest.tsv \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path ${prodir}/outputs/qiime2/qza/all.qza \
--type 'SampleData[PairedEndSequencesWithQuality]'

### Summarize read depth and quality scores
qiime demux summarize \
--i-data ${prodir}/outputs/qiime2/qza/all.qza \
--o-visualization ${prodir}/outputs/qiime2/qzv/all.qzv

#bash ${prodir}/bash/qiime2/all_dada2_qiime2.sh
