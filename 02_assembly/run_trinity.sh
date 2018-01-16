#!/bin/bash

bamFile=mapping.bam
maxIntron=10000
maxMemory=150G
maxCPU=32

/storage/software/intel64/trinityrnaseq_v2.2.0/Trinity \
    --genome_guided_bam "$bamFile" \
    --genome_guided_max_intron "$maxIntron" \
    --max_memory "$maxMemory" \
    --CPU "$maxCPU" 
