from utils import parse_fasta
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd


if __name__ == "__main__":

    model_path = Path("/storage/group/izg5139/default/Avoidmers-DNA/models")
    # model_path = Path("model_organisms")

    def extract_name(f: str | Path) -> str:
        f = Path(f).name
        if f.startswith("GC"):
            return '_'.join(f.split("_")[:2])
        return f.split(".fa")[0]

    model_organisms = [f for f in model_path.glob("*.fa")]
    model_organisms = model_organisms + [f for f in model_path.glob("*.fna")]
    data = defaultdict(list)
    original_sizes = defaultdict(list)
    genome_size = {}

    for model in model_organisms:
        model_name = extract_name(model)
        print(f"Processing model organism {model}...")
        total = 0
        for seqID, seq in parse_fasta(model):
            sequence_length = len(seq)
            data["organism"].append(model_name)
            data["seqID"].append(seqID)
            data["length"].append(sequence_length)
            total += sequence_length

        genome_size.update({model_name: total})
        chrom_size = model.parent.parent.joinpath('model_organsisms', model_name + 'chrom.sizes')

        if chrom_size.is_file():
            with chrom_size.open('r') as f:
                for line in f:
                    line = line.strip().split("\t")
                    seqID, size = line
                    original_sizes["organism"].append(model_name)
                    original_sizes["seqID"].append(seqID)
                    original_sizes["lengthOri"].append(int(size))

    with open("genome_sizes.txt", mode="w") as f:
        for model, total in genome_size.items():
            f.write(f"{model}\t{total}\n")

    df = pd.DataFrame(data)
    original_sizes = pd.DataFrame(original_sizes)
    if chrom_size.is_file():
        df_merged = df.merge(original_sizes, on=["organism", "seqID"], how="left")
        df_merged.loc[:, "diff"] = np.abs(df_merged["length"] - df_merged["lengthOri"])
        assert df_merged[df_merged['diff'] > 0].shape[0] == 0

        df_merged.to_csv("genome_chromosome_sizes.txt", sep="\t", mode="w", index=False)





