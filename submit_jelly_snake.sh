#!/bin/bash

j=${1:-1}

if [[ ! -n ${SSH_CONNECTION} ]];
then
    snakemake --snakefile scripts/jelly_snake.smk --cores 1 --latency-wait 5 --keep-going
else
    snakemake --use-conda --snakefile scripts/jelly_snake.smk \
	    --rerun-incomplete \
	    --keep-incomplete \
	    --rerun-triggers mtime \
	    --keep-going \
	    --cores $j \
	    --latency-wait 45 \
	    --cluster-config config/cluster_settings.yaml \
	    --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOut/{cluster.jobName}-%j.out -J {cluster.jobName} -e jobOut/{cluster.jobName}-%j.err" -j $j
fi


