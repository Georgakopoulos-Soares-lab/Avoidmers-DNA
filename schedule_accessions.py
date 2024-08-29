from pathlib import Path
from utils import extract_name, parse_fasta
import argparse
from collections import defaultdict
import pandas as pd


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--accession_path", type=str, default="accessions")

    args = parser.parse_args()
    accession_path = Path(args.accession_path).resolve()

    accession_mapper = defaultdict(list)

    for accession in accession_path.glob("*.gz"):

        accession_name = extract_name(accession)

        for seqID, _ in parse_fasta(accession):

            accession_mapper["accession"].append(str(accession))
            accession_mapper["accession_name"].append(accession_name)
            accession_mapper["seqID"].append(seqID)


    accession_mapper = pd.DataFrame(accession_mapper)
    accession_mapper.to_csv("accession_seqID_mapping.txt", mode="w", sep="\t", index=False, header=True)






