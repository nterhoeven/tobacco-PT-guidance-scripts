#!/bin/bash

SeqChunker --chunk-number 10 --out proteins-chunk.%02d.fa proteins_combined-IUPAC-len10.fa
