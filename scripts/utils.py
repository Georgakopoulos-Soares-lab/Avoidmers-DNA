import gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser
from pathlib import Path
import os
from typing import Iterable

def parse_fasta(accession: str | os.PathLike[str]) -> Iterable[tuple[str, str]]:
    accession = Path(accession)
    if accession.name.endswith(".gz"):
        file = gzip.open(accession, 'rt')
    else:
        file = open(accession, 'r')

    for seqID, seq in SimpleFastaParser(file):
        yield seqID.split(' ')[0], seq.replace("\n", "").lower()

    file.close()

def search_chromosome(accession: str | os.PathLike[str], target_seqID: str) -> tuple[str, str]:
    for seqID, seq in parse_fasta(accession):
        if target_seqID == seqID:
            return seqID, seq
    raise ValueError(f"Failed to detect chromosome {target_seqID}.")

def extract_name(accession: str) -> str:
    accession = Path(accession).name
    if ".gz" in accession:
        accession = accession.split(".gz")[0]
    return ".".join(accession.split(".")[:-1])

def merged_seq(row: dict) -> str:
    starts = list(map(int, row["allStarts"].split("|")))
    sequences = row["allSequences"].split("|")
    assert int(row["overlapCount"]) == len(sequences) == len(starts), f"Invalid sequence count for {row} during merge operation."
    merged_sequence = ""
    if len(starts) == 1:
        return sequences[0]

    for i, (s1, s2) in enumerate(zip(starts, starts[1:])):
        merged_sequence += sequences[i][:s2-s1]

    merged_sequence += sequences[-1]
    assert len(merged_sequence) == int(row["end"]) - int(row["start"]), f"Invalid sequence merge {row}"
    return merged_sequence


def reverse_complement(n: str) -> str:
    match n:
        case 'a':
            return 't'
        case 't':
            return 'a'
        case 'c':
            return 'g'
        case 'g':
            return 'c'
        case _:
            raise ValueError(f'Invalide nucleic acide {n}.')


def reverse_kmer(kmer: str) -> str:
    return ''.join(reverse_complement(n) for n in kmer)[::-1]
