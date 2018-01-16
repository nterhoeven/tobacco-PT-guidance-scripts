#!/bin/bash


threads=32

for i in *R1_001.trimmed.fastq.gz
do 
    first="$i"
    second=$(echo "$i" | sed -e 's/R1/R2/')
    base=$(basename "$first" _R1_001.trimmed.fastq.gz)

    echo "###########################################"
    echo "running sample $base"
    echo "$first -> $second -> $base"

    /software/Salmon-0.7.2/bin/salmon quant -i transcriptome-index -l A -1 "$first" -2 "$second" -p "$threads" -o quants/"$base"_quant
done