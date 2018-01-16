#!/bin/bash

threads=32

for i in Nt*_R1*
do
    echo "####################################################"
    date
    echo "starting trimming for $i"

    reads="$i"
    mates=$(echo "$i" | sed -e 's/R1/R2/')
    
    readsOut=$(echo "$reads" | sed -e 's/\.fastq\.gz/\.trimmed\.fastq\.gz/')
    readsOutSingle=$(echo "$reads" | sed -e 's/\.fastq\.gz/\.trimmed\.single\.fastq\.gz/')
    matesOut=$(echo "$mates" | sed -e 's/\.fastq\.gz/\.trimmed\.fastq\.gz/')
    matesOutSingle=$(echo "$mates" | sed -e 's/\.fastq\.gz/\.trimmed\.single\.fastq\.gz/')

    echo "using filenames:"
    echo "reads: $reads"
    echo "mates: $mates"
    echo "reads out: $readsOut"
    echo "reads out single: $readsOutSingle"
    echo "mates out: $matesOut"
    echo "mates out single: $matesOutSingle"

    echo "starting trimmomatic"
    java -jar /software/trimmomatic/trimmomatic-0.33.jar PE -threads "$threads" -phred33 "$reads" "$mates" "$readsOut" "$readsOutSingle" "$matesOut" "$matesOutSingle" ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:95

    echo "done"
done

echo "finished all"