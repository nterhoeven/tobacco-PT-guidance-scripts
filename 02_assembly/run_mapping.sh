#!/bin/bash

export PATH="/software/bowtie2/:$PATH"
reads1=$(echo *R1_001.trimmed.fastq.gz | sed -e 's/ /,/g')
reads2=$(echo *R2_001.trimmed.fastq.gz | sed -e 's/ /,/g')
genome=Ntab-TN90_AYMY-SS
threads=32

#run tophat
echo "############################"
date
echo "running tophat"
echo "############################"

/software/tophat/tophat -p "$threads" "$genome" "$reads1" "$reads2"

echo "############################"
date
echo "finished tophat"
echo "############################"

