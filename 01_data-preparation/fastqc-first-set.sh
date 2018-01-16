#!/bin/bash

for i in Nt*.fastq.gz
do
    echo "running fastqc for $i"
    /software/fastqc/fastqc --threads 32 "$i"
done
