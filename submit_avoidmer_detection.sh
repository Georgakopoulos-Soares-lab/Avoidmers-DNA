#!/bin/bash

#SBATCH --mem=20GB
#SBATCH --time=48:00:00
#SBATCH --job-name=AvoidmerDetection
#SBATCH --output=logs/%x_%j.err
#SBATCH --error=logs/%x_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --partition=open


protocol=${1:-"abacaba"}
kmer_length=${2:-13}
alphabet=${3:-"agct"}
save_kmers=${4:-1}
outdir=${5:-"expected_avoidmers"}

python generate_avoidmers.py --kmer_length ${kmer_length} --outdir ${outdir} --alphabet ${alphabet} --save_kmers ${save_kmers} --protocol ${protocol}
