#!/bin/bash

for i in proteins-chunk.*.fa
do
    sbatch -c 30 --mem 100G ./run_Interproscan.sh "$i"
done

    
