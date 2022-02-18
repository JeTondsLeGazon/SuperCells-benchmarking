#!/bin/bash

echo "Dowloading data"
mkdir -p data/

curl -s https://zenodo.org/record/5048449/files/bulk_rnaseq.tar.gz?download=1 --output data/bulk_rnaseq.tar.gz
curl -s https://zenodo.org/record/5048449/files/sc_rnaseq.tar.gz?download=1 --output data/sc_rnaseq.tar.gz

echo "Successfully downloaded data files"

tar -xf data/bulk_rnaseq.tar.gz -C data/
tar -xf data/sc_rnaseq.tar.gz -C data/

rm data/bulk_rnaseq.tar.gz
rm data/sc_rnaseq.tar.gz

echo "Data available in data folder"
