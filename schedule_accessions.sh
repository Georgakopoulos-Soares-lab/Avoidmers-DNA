#!/bin/bash

dir=${1:-"model_organisms_local"}
accessions=`(find $(pwd)/${dir}/ -maxdepth 1 -type f ! -name "*.gff")`

mkdir -p scheduling/
echo -e accession'\t'accession_name'\t'seqID >> scheduling/scheduled_accessions.txt
for accession in ${accessions[@]};
do
    accession_name=$(basename "${accession%.*}")
    seqkit fx2tab -n ${accession} > scheduling/${accession_name}.seqID.txt
    awk -v prefix=${accession_name} -v name=${accession} '{ print name "\t" prefix "\t" $1 }' scheduling/${accession_name}.seqID.txt >> scheduling/scheduled_accessions.txt
done

mv scheduling/scheduled_accessions.txt .
rm -rf scheduling
