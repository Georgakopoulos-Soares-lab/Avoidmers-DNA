#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --mem=16GB

# source activate my_env

bash submit_zimin_snake.sh 100
