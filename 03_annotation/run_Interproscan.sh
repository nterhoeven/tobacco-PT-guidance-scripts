#!/bin/bash

echo "##################################"
date
echo "starting Interproscan for chunk $1"
echo "##################################"

/storage/compevolbiol/software/Interproscan/interproscan-5.25-64.0/interproscan.sh --cpu 30 --goterms --input $1
