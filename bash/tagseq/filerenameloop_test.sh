#!/bin/bash
testletters=(a b c d e f)
testletternumbers=(a1 b2 c3 d4 d5 d6)

for i in {0..5}
do
touch ${testletters[i]}.txt
done

for i in {0..5}
do
mv -- ${testletters[i]}.txt ${testletternumbers[i]}.txt
done

readfiles=(ls ${prodir}/data/quantseq_reads)
quantseq_samples=(cat quantseq_samples.txt)

for i in {0..60}
do
mv -- ${readfiles[i]} ${quantseq_samples[i]}.fastq.gz
done

for file in *.nii
do read line
   mv -v "${file}" "${line}"
done < list
