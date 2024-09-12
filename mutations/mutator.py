import re
import abc
from pathlib import Path
import csv

class Mutator(abc.ABC):

    def __init__(self, alphabet: str = "agct") -> None:
        self.alphabet = alphabet
        # self.avoidmers = set()

    @abc.abstractmethod
    def is_invariant(self, sequence: str) -> bool:
        raise NotImplementedError()

    def mutate(self, sequence: str) -> dict[str, float | str | list[int]]:
        snp_invariance, snp_positions = self.mut_snp(sequence)
        del_invariance, del_positions = self.mut_deletions(sequence)
        ins_invariance, ins_positions = self.mut_insertions(sequence)

        return {
                "snp_invariance": snp_invariance,
                "snp_positions": snp_positions,
                "del_invariance": del_invariance,
                "del_positions": del_positions,
                "ins_invariance": ins_invariance,
                "ins_positions": ins_positions,
                "sequence": sequence,
                "is_invariant": self.is_invariant(sequence),
            }

    def mut_snp(self, sequence: str | list) -> tuple[float, list[int]]:
        invariant = 0
        N = len(sequence)
        total_mutations = 3 * N
        sequence = list(sequence)
        positions = [0 for _ in range(N)]

        # 3 X N total SNP
        for i in range(N):
            allele = sequence[i] # allele

            for mut in self.alphabet: # agct = nucleotide alphabet
                if mut == allele:  # continue
                    continue

                sequence[i] = mut
                mutated_sequence = ''.join(sequence)

                # caching
                # if mutated_sequence in self.avoidmers:
                #     positions[i] += 1
                #     survived += 1
                #     total += 1
                #     sequence[i] = allele
                #     continue

                # bottleneck
                if self.is_invariant(mutated_sequence):
                    invariant += 1
                    positions[i] += 1
                    # self.avoidmers.add(seq)

            sequence[i] = allele

        return invariant / float(total_mutations), positions

    def mut_insertions(self, sequence: str) -> tuple[float, list[int]]:
        invariant = 0
        N = len(sequence)
        total_mutations = (N+1) *   4
        positions = [0 for _ in range(N+1)]

        for i in range(N+1):
            for mut in self.alphabet: # agct = nucleotides
                mutated_sequence = sequence[:i] + mut + sequence[i:]

                # caching
                # if mutated_sequence in self.avoidmers:
                #     survived += 1
                #     total += 1
                #     position[i] += 1
                #     continue

                # bottleneck
                if self.is_invariant(mutated_sequence):
                    invariant += 1
                    positions[i] += 1
                    # self.avoidmers.add(mutated_sequence)

        return invariant / float(total_mutations), positions

    def mut_deletions(self, sequence: str) -> tuple[float, list[int]]:
        invariant = 0
        N = len(sequence)
        total_mutations = N
        positions = [0 for _ in range(N)]

        # N deletions
        for i in range(N):
            mutated_sequence = sequence[:i] + sequence[i+1:]

            # caching
            # if mutated_sequence in self.avoidmers:
            #     survived += 1
            #     total += 1
            #     positions[i] += 1
            #     continue

            # bottleneck
            if self.is_invariant(mutated_sequence):
                positions[i] += 1
                invariant += 1
                # self.avoidmers.add(mutated_sequence)

        return invariant / float(total_mutations), positions

    def __repr__(self) -> str:
        return f"Sequence Mutator: {type(self).__name__} using Alphabet({''.join(self.alphabet)})."

class AvoidmerMutator(Mutator):

    def __init__(self, alphabet: str = "agct") -> None:
        super().__init__(alphabet)

    def is_invariant(self, sequence: str) -> bool:
        return re.search(r"(.+)(.+)\1(.+)\1\2\1", sequence) is None

class AbaMutator(Mutator):

    def is_invariant(self, sequence: str) -> bool:
        return re.search(r"(.+)(.+)\1", sequence) is None


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", type=str, default="sequences.txt")
    parser.add_argument("--out", type=str, default="sequence_mutations.mi")
    parser.add_argument("--alphabet", type=str, default="agct")

    args = parser.parse_args()
    sequences = Path(args.sequences).resolve()
    out = Path(args.out)
    alphabet = args.alphabet

    if not sequences.is_file():
        raise FileNotFoundError(f"File {sequences} does not exist.")

    mutator = AvoidmerMutator(alphabet=alphabet)
    MI_FIELDS = ["sequence", "is_invariant", "snp_invariance", "snp_positions", "del_invariance", "del_positions", "ins_invariance", "ins_positions"]
    mutation_types = ["snp", "del", "ins"]
    with out.open("w") as g, sequences.open("r") as f:
        writer = csv.DictWriter(g, fieldnames=MI_FIELDS, delimiter="\t")
        for sequence in f:
            sequence = sequence.strip()
            invariance_data = mutator.mutate(sequence)
            for pos in mutation_types:
                invariance_data[f"{pos}_positions"] = ','.join(str(p) for p in invariance_data[f"{pos}_positions"])

            writer.writerow(invariance_data)
