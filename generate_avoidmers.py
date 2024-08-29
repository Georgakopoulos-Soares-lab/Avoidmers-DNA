from detection import get_search_protocol
from pathlib import Path

class Avoidmers:

    def __init__(self, protocol: str, alphabet: str = "agct") -> None:
        self.alphabet = alphabet
        self.alphabet_size = len(alphabet)
        self.protocol = protocol
        self.is_free = get_search_protocol(protocol)
        self.total_avoidmers = {1: self.alphabet_size}
        self.probabilities = {1: 1}

    def generate_avoidmers(self, kmer_length: int) -> list[str]:
        if kmer_length == 1:
            return list(self.alphabet)

        avoidmers = self.generate_avoidmers(kmer_length-1)
        kmers = [avoidmer + n for avoidmer in avoidmers for n in self.alphabet if self.is_free(avoidmer+n)]
        self.total_avoidmers.update({kmer_length: len(kmers)})
        self.probabilities.update({kmer_length: (1e2 * self.total_avoidmers[kmer_length]) / self.alphabet_size ** kmer_length})
        return kmers

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="""Avoidmer Detector.""")
    parser.add_argument("--kmer_length", type=int, default=8)
    parser.add_argument("--alphabet", type=str, default="agct")
    parser.add_argument("--outdir", type=str, default="avoidmers_expected_vs_observed")
    parser.add_argument("--protocol", type=str, default="abacaba", choices=["square", "aba", "abacaba"])
    parser.add_argument("--save_kmers", type=int, choices=[0, 1], default=1)

    args = parser.parse_args()
    kmer_length = args.kmer_length
    alphabet = args.alphabet
    protocol = args.protocol
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(exist_ok=True)

    save_kmers = args.save_kmers

    avoidmer = Avoidmers(alphabet=alphabet, protocol=protocol)
    if save_kmers:
        kmers = avoidmer.generate_avoidmers(kmer_length=kmer_length)
        with open(f"{outdir}/avoidmers_{kmer_length}_{alphabet}_{protocol}.txt", mode="w") as f:
            for kmer in kmers:
                f.write(f"{kmer}\n")

    total_avoidmers = avoidmer.total_avoidmers
    probabilities = avoidmer.probabilities

    with open(f"{outdir}/expected_probabilities_{kmer_length}_{alphabet}_{protocol}.txt", mode="w") as f:
        for kmer_length in probabilities:
            f.write(f"{kmer_length}\t{len(alphabet) ** kmer_length}\t{total_avoidmers[kmer_length]}\t{probabilities[kmer_length]}\n")
