import numpy as np
from dataclasses import dataclass, field
from pathlib import Path
import json

@dataclass
class Scheduler:

    suffix: str = field(default=".txt")

    def schedule(self, indir: str | Path, total_buckets: int = 100, fout: str | Path = "schedule.json") -> str | Path:
        indir = Path(indir).resolve()
        if indir.is_dir():
            items = [str(f) for f in indir.glob(f".{self.suffix}")]
        elif indir.is_file():
            items = []
            with open(indir, 'r') as f:
                for item in f:
                    item = item.strip()
                    items.append(item)
        elif iter(indir):
            items = indir

        jobs = {i: job.tolist() for i, job in enumerate(np.array_split(items, total_buckets))}
        with open(fout, mode='w') as f:
            json.dump(jobs, f, indent=4)

        return fout


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--suffix", type=str, default=".txt")
    parser.add_argument("--indir", type=str)
    parser.add_argument("--total_buckets", type=int, default=100)
    parser.add_argument("--fout", type=str, default="schedule.json")

    args = parser.parse_args()
    suffix = args.suffix
    indir = args.indir
    fout = args.fout
    total_buckets = args.total_buckets

    scheduler = Scheduler(suffix=suffix)
    scheduler.schedule(indir=indir,
                       total_buckets=total_buckets,
                       fout=fout
                       )
