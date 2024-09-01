from detection import get_search_protocol
from generate_avoidmers import Avoidmers
from pathlib import Path
import subprocess
import platform
import os
from termcolor import colored
from validate_counts import count_kmers

if "SSH_CONNECTION" in os.environ:
    print("Initializing snakemake. Using remote configuration...")
    configfile: 'config/config_jelly_server.yaml'
else:
    print("Initializing snakemake. Using local configuration...")
    configfile: 'config/config_jelly.yaml'


indir = config['indir']
outdir = Path(config['outdir']).resolve()
outdir.mkdir(exist_ok=True)
protocol = config['protocol']
max_kmer_length = int(config['max_kmer_length'])

print(f"INDIR ---------> {indir}")
print(f"OUTDIR --------> {outdir}")
print(f"PROTOCOL -------> {protocol}")
print(f"MAX KMER LENGTH ----> {max_kmer_length}")

is_free = get_search_protocol(protocol)
ALL = glob_wildcards('%s/{accession}.{suffix}' % indir)
SAMPLES = ALL.accession
SUFFIXES = ALL.suffix

print(f"TOTAL SAMPLES DETECTED: {len(SAMPLES)}.")
print("INITIALIZING SNAKEMAKE...")

sample_suffix_mapping = {sample: suffix for sample, suffix in zip(SAMPLES, SUFFIXES)}

rule all:
    input:
        expand('%s/{accession}.avoidmer_probability.txt' % outdir, accession=SAMPLES)

rule count_kmers:
    input:
        lambda wc: '%s/%s.%s' % (indir, wc.accession, sample_suffix_mapping[wc.accession])
    output:
        '%s/{accession}.fasta_{kmer_length}mers' % outdir
    params:
        validate_counts=int(config['validate_counts']),
        genome_sizes=config['genome_sizes']
    run:
       suffix = sample_suffix_mapping[wildcards.accession]
       infile = input[0]
       outfile = f"{outdir}/{wildcards.accession}.raw.fasta_{wildcards.kmer_length}mers"
       out_dumped = f"{outdir}/{wildcards.accession}.fasta_{wildcards.kmer_length}mers"

       # HONOR CONTRACT
       assert wildcards.kmer_length.isdigit(), f"Invalid kmer length {wildcards.kmer_length}"
       kmer_length = int(wildcards.kmer_length)
       subprocess.run(f"jellyfish count -m {kmer_length} -s 100M -t 10 -o {outfile} {infile}", check=True, shell=True)
       subprocess.run(f"jellyfish dump -c {outfile} > {out_dumped}", check=True, shell=True)
       os.remove(outfile)

       genome_sizes = {}
       with open(params.genome_sizes, 'r') as f:
           for line in f:
               line = line.strip()
               accession, total = line.split("\t")
               total = int(total)
               genome_sizes.update({accession: total})

       if params.validate_counts and int(wildcards.kmer_length) < 3:
           kmer_python_counts = count_kmers(infile, int(wildcards.kmer_length))
           f = open(f"{outdir}/{wildcards.accession}.python.fasta_{wildcards.kmer_length}mers", mode="w")
           for kmer, occurrences in kmer_python_counts.items():
               f.write(f"{kmer} {occurrences}\n")
           f.close()
           del kmer_python_counts['total_skipped']
           
           jelly_counts = {}
           total_occurrences = 0 
           with open(out_dumped, 'r') as f:
               for line in f:
                   line = line.strip().split(" ")
                   kmer, occurrences = line
                   kmer = kmer.lower()
                   occurrences = int(occurrences)
                   jelly_counts.update({kmer: occurrences})
                   total_occurrences += occurrences
           if wildcards.accession in genome_sizes:
               total = genome_sizes[wildcards.accession]
               total_kmers  = total - int(wildcards.kmer_length) + 1
               print(f'{input[0]}> CHECKPOINT: >GENOME SIZE IS CORRECT< VALIDATION STATUS: {total},{total_occurrences},{total_kmers}.')
          

           try:
               assert set(kmer_python_counts.keys()) == set(jelly_counts.keys())
               print(colored(f"{input[0]}> CHECKPOINT: >KMER SETS ARE EQUAL< VALIDATION STATUS: SUCCEEDED.", "green"))
           except AssertionError:
               print(colored(f"{input[0]}> CHECKPOINT: >KMER SETS ARE EQUAL< VALIDATION STATUS: FAILED.", "red"))

           try:
               for kmer in kmer_python_counts:
                   assert kmer_python_counts[kmer] == jelly_counts.get(kmer, 0)
               print(colored(f"{input[0]}> CHECKPOINT: >KMER COUNTS ARE THE SAME< VALIDATION STATUS: SUCCEEDED.", "green"))
           except AssertionError:
               print(colored(f"{input[0]}> CHECKPOINT: >KMER COUNTS ARE THE SAME< VALIDATION STATUS: FAILED.", "red"))



rule find_avoidmers:
    input:
        expand('%s/{{accession}}.fasta_{kmer_length}mers' % outdir, kmer_length=range(1, max_kmer_length))
    output:
        '%s/{accession}.avoidmer_probability.txt' % outdir
    run:
        avoidmer_probabilities = {}
        avoidmer_occurrences = {}
        kmer_set = {}
        for kmer_length in range(1, max_kmer_length):
            infile = f'{outdir}/{wildcards.accession}.fasta_{kmer_length}mers'
            total_kmers = 0
            total_avoidmers = 0
            with open(infile, 'r') as f:
                for line in f:
                    line = line.strip().split(" ")
                    kmer, occurrences = line
                    kmer = kmer.lower()
                    occurrences = int(occurrences)

                    is_avoidmer = is_free(kmer)
                    total_kmers += occurrences

                    if is_avoidmer:
                        total_avoidmers += occurrences

            probability_class_1 = (1e2 * total_avoidmers) / (total_kmers)
            kmer_set.update({kmer_length: total_kmers})
            avoidmer_probabilities.update({kmer_length: probability_class_1})
            avoidmer_occurrences.update({kmer_length: total_avoidmers})

        with open(output[0], mode='w') as f:
            for kmer_length, probability in avoidmer_probabilities.items():
                f.write(f"{kmer_length}\t{kmer_set[kmer_length]}\t{avoidmer_occurrences[kmer_length]}\t{probability}\n")
