#!/bin/bash

mkdir -p model_organisms
cd model_organisms/

fin=${1:-model_organisms}
while read model;
do
  echo "Downloading model organism ${model}."
  if [[ ! -f ${model}.fa.gz ]];
  then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/${model}/bigZips/${model}.fa.gz
  fi

  if [[ ! -f ${model}.chrom.sizes ]];
  then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/${model}/bigZips/${model}.chrom.sizes
  fi

  if [[ ! -f ${model}.md5sum.txt ]];
  then
    wget -O ${model}.md5sum.txt https://hgdownload.soe.ucsc.edu/goldenPath/${model}/bigZips/md5sum.txt
  fi

  if [[ ! -f ${model}.ncbiRefSeq.gtf.gz ]];
  then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/${model}/bigZips/genes/${model}.ncbiRefSeq.gtf.gz
  fi
done < "../$fin"

echo "Downloading remaining NCBI genomes..."
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gff.gz

