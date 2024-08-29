from utils import parse_fasta, reverse_kmer
from pathlib import Path
from collections import defaultdict

def count_kmers(accession, kmer_length):
    counts = defaultdict(int)
    nucleotides = {"a", "g", "c", "t"}
    for _, seq in parse_fasta(accession):
        for i in range(len(seq)-kmer_length+1):
            chunk = seq[i:i+kmer_length]

            if any(c not in nucleotides for c in chunk):
                continue
            counts[chunk] += 1

    return dict(counts)

def extract_name(f):
    f = Path(f).name
    if f.startswith("GC"):
        return '_'.join(f.split("_")[:2])
    return f.split(".fa")[0] if ".fa" in f else f.split(".fna")[0]

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--accession", type=str)
    parser.add_argument("--kmer_length", type=int, default=3)
    parser.add_argument("--jellycounts", type=str)

    args = parser.parse_args()
    accession = args.accession
    kmer_length = args.kmer_length

    jellydir = Path("jelly")

    # read jelly counts
    name = extract_name(accession)
    jellypath = jellydir.joinpath(name, name + f".unzipped.fasta_{kmer_length}mers.dumped")
    jelly_counts = {}
    total_kmers = 0

    print(f"Total kmers detected {total_kmers} out of {4 ** kmer_length}.")
    counts = count_kmers(accession, kmer_length)
    breakpoint()

    with open(jellypath, mode='r') as f:
        for line in f:
            line = line.strip().split(" ")
            kmer, occurrences = line
            kmer = kmer.lower()
            jelly_counts.update({kmer: int(occurrences)})
            reverse_com = reverse_kmer(kmer)
            jelly_counts.update({reverse_com: int(occurrences)})
            total_kmers += len({reverse_com, kmer})


    assert set(counts.keys()) == set(jelly_counts.keys())
    for key in counts:
        assert counts[key] == jelly_counts[key], f"Error at {key}! {counts[key]} != {jelly_counts[key]}."

    print("ALL GOOD!")
