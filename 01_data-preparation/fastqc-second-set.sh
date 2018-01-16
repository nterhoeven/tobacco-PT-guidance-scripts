#!/bin/bash

for i in *.fastq.gz
do
    echo "running fastqc for $i"
    outdir="fastqc_$i"
    mkdir "$outdir"
    /software/fastqc/fastqc --threads 32 --outdir "$outdir" "$i"
done
