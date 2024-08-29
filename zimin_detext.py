from utils import search_chromosome, extract_name
import os
import csv
from pathlib import Path
from detection import get_search_protocol
import pandas as pd

def detect_pattern(accession: os.PathLike[str],
                            seqID: str,
                            kmer_length: int,
                            search_protocol: str,
                            out: os.PathLike[str],
                            ) -> None:

    out = Path(out).resolve()
    out.mkdir(exist_ok=True)
    outdir = out.joinpath(f"pattern_extractions_{search_protocol}")
    outdir.mkdir(exist_ok=True)

    is_pattern_free = get_search_protocol(search_protocol)
    kmer_length = int(kmer_length)
    seqID, seq = search_chromosome(accession, seqID)
    seq_length = len(seq)

    print(f"Initializing detection algorithm with search protocol {search_protocol} for accession {accession} on chromosome {seqID} for k-mer length {kmer_length}...", flush=True)
    print(f"Total sequence length {seqID}:{seq_length}.\nRedirecting output to {outdir}...", flush=True)

    accession_name = extract_name(accession)
    filename = f"{accession_name}_{search_protocol}_words_length_{kmer_length}_seq_{seqID}.txt"

    out = open(outdir.joinpath(filename), mode="w")
    writer = csv.DictWriter(out,
                            fieldnames=["seqID", "start", "end", "sequence", "length", "type"],
                            delimiter=","
                            )

    seq_length = len(seq)
    nucleotides = {"a", "g", "c", "t"}
    from tqdm import tqdm

    for i in tqdm(range(seq_length-kmer_length+1), total=seq_length-kmer_length+1):
        for lookahead in range(kmer_length, seq_length-i+1):
            chunk = seq[i:i+lookahead]

            # i = 0 | lookahead = 8
            # AGAGAGA|U]X
            # i = 1 | lookahead = 7
            #
            #  GAGAGA|U]X

            if any(c not in nucleotides for c in chunk):
                # break since it contains unsequenced data
                break

            # check if chunk is avoidmer
            is_avoidmer = is_pattern_free(chunk)

            # if it is not, break and move on 1bp to new chunk
            if lookahead == kmer_length and not is_avoidmer:
                break
            # if lookahead is higher than kmer-length, it means that we are attempting extension
            # if it's not avoidmer anymore, go back 1bp and register as the maximal avoidmer
            elif not is_avoidmer:
                info = {
                        "seqID": seqID,
                        "start": i,
                        "end": i+lookahead-1,
                        "sequence": chunk[:-1],
                        "length": len(chunk[:-1]),
                        "type": search_protocol
                    }

                writer.writerow(info)
                assert lookahead-1 == len(chunk[:-1])
                # break since current k-mer is no longer extendable
                break
        else:
            # search continued at the end of the sequence without breaking
            # thus, k-mer was maximally extendable avoidmer
            info = {
                    "seqID": seqID,
                    "start": i,
                    "end": seq_length,
                    "sequence": chunk,
                    "length": len(chunk),
                    "type": search_protocol
                    }
            writer.writerow(info)
    out.close()

    # A|GCAGA|TATAGGACTA(U)
    #   GCAGA|TATAGGACTA(U)
    #    CAGA|TATAGGACTA

    avoidmers_df = pd.read_csv(outdir.joinpath(filename))\
                    .groupby(["seqID", "end"], as_index=False)\
                    .agg({
                          "start": "min",
                          "sequence": "first",
                          "length": "max",
                          "type": "first"
                          })

    maximal_filename = f"{accession_name}_{search_protocol}_words_length_{kmer_length}_seq_{seqID}.maximal.txt"
    avoidmers_df.to_csv(outdir.joinpath(maximal_filename), sep=",", mode="w", index=False)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="""Pattern Avoidance Detector|Fixed k-mer length""")
    parser.add_argument("--accession", type=str)
    parser.add_argument("--seqID", type=str, default="chr1")
    parser.add_argument("--kmer_length", type=int, default=60)
    parser.add_argument("--out", type=str, default="patterns")
    parser.add_argument("--search_protocol",
                        type=str,
                        default="abacaba",
                        choices=["aba", "abacaba", "square"])

    args = parser.parse_args()
    accession = args.accession
    seqID = args.seqID
    kmer_length = args.kmer_length
    out = args.out
    search_protocol = args.search_protocol

    detect_pattern(
                accession=accession,
                seqID=seqID,
                kmer_length=kmer_length,
                search_protocol=search_protocol,
                out=out
            )
