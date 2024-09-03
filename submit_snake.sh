#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=AvoidmerSnake
#SBATCH --out=jobOut/%j_%x.out
#SBATCH --error=jobOut/%j_x.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

# source activate my_env

jobs=${1:-100}
echo $jobs
bash submit_zimin_snake.sh $jobs
