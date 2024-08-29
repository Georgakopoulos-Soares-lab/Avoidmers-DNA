#!/bin/bash


protocols = ("aba" "abacaba" "square")

for protocol in ${protocols[@]};
do
    sbatch submit_avoidmer_detection.sh ${protocol}
done
