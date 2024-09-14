import numpy as np
import pandas as pd
from dataclasses import dataclass

@dataclass
class WindowMaker:

    window_size: int = 500
    base: int = 1

    def make_windows(self, df: pd.DataFrame, loci: str) -> pd.DataFrame:

        if loci != "start" and loci != "end":
            raise ValueError(f"Unknown loci {loci}.")

        if "chromosome" in df:
            df = df.rename(columns={"chromosome": "seqID"})

        if "seqID" not in df:
            raise KeyError("seqID column is not present in the dataframe.")

        if self.base == 1:
            df[loci] = df[loci] - 1

        if loci == "start":
            df["end"] = df[loci] + self.window_size + 1
            df["start"] = np.maximum(df[loci] - self.window_size, 0)
        elif loci == "end":
            df["start"] = np.maximum(df[loci] - self.window_size, 0)
            df["end"] = df[loci] + self.window_size + 1

        columns = [col for col in df.columns if col != "start" and col != "end" and col != "seqID"]
        return df[["seqID", "start", "end"] + columns]
