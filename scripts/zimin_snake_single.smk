from pathlib import Path
import pybedtools
from pybedtools import BedTool
from zimin_detext import detect_pattern
from utils import extract_name, merged_seq
import pandas as pd

if "SSH_CONNECTION" in os.environ:
    print("Initializing snakemake. Using remote configuration...")
    configfile: 'config/config_server.yaml'
else:
    print("Initializing snakemake. Using local configuration...")
    configfile: 'config/config.yaml'

out = Path(config['out']).resolve()
protocol = config['protocol']
kmer_length = int(config['kmer_length'])
# samples, = glob_wildcards("%s/{assembly}.fa.gz" % data_path)

info = pd.read_table(config['info'], skiprows=1, names=["accession", "name", "seqID"], header=None)
info_assembly = info.set_index("name")["accession"].to_dict()

rule all:
    input:
        expand('%s/pattern_merged_%s/{assembly}_%s_words_length_%s.merged.txt' % (out, protocol, protocol, kmer_length),
               assembly=info['name'].drop_duplicates()
               )

rule extractPatterns:
    input:
        assembly_location=lambda wc: info_assembly[wc.assembly]
    output:
        '%s/pattern_extractions_%s/{assembly}_%s_words_length_%s_seq_{seqID}.txt' % (out, protocol, protocol, kmer_length)
    params:
        kmer_length=int(config['kmer_length']),
        protocol=config['protocol']
    run:
        detect_pattern(accession=input[0],
                       seqID=wildcards.seqID,
                       kmer_length=params.kmer_length,
                       search_protocol=params.protocol,
                       out=out
                    )

rule mergeExtractions:
    input:
        lambda wc: expand('%s/pattern_extractions_%s/{{assembly}}_%s_words_length_%s_seq_{seqID}.txt' % (out, protocol, protocol, kmer_length),
                seqID=info[info['name'] == wc.assembly]['seqID']
               )
    output:
        '%s/pattern_merged_%s/{assembly}_%s_words_length_%s.merged.txt' % (out, protocol, protocol, kmer_length)
    params:
        out=Path(config['out']).resolve(),
        protocol=config['protocol'],
        kmer_length=int(config['kmer_length'])
    run:
        Path(config['temp_dir']).mkdir(exist_ok=True)
        pybedtools.helpers.set_tempdir(config['temp_dir'])
        pybedtools.set_bedtools_path(config['bedtools'])

        df = []
        # '%s/pattern_merged_%s/{assembly}_%s_words_length_%s.extended.merged.txt' % (out, protocol, protocol, kmer_length)
        for extracted_seq in input:
            print(extracted_seq)
            temp_df = pd.read_csv(extracted_seq)
            df.append(temp_df)

        df = pd.concat(df, axis=0)
        df.to_csv(output[0], sep=",", index=False, mode="w")

        # save maximum Zimin-Avoiding Sequences
        filename = f"{params.out}/pattern_merged_{params.protocol}/{wildcards.assembly}_{params.protocol}_words_length_{params.kmer_length}.maximum_zimin.txt"
        df.groupby("seqID").agg({"length": "max"})\
                            .to_csv(
                                filename,
                                mode="w",
                                index=True,
                                sep=","
                            )

        print(df)
        # save merged coordinates
        df_merged = pd.read_table(
                   BedTool.from_dataframe(df)\
                             .sort()\
                             .merge(c=["2", "4", "4"], o=["collapse", "collapse", "count"], delim="|")\
                             .fn,
                     header=None,
                     names=["seqID", "start", "end", "allStarts", "allSequences", "overlapCount"]
                    )

        df_merged.loc[:, "mergedSequence"] = df_merged[["allStarts", "allSequences"]].apply(lambda row: merged_seq(row), axis=1)
        df_merged.drop(columns=["allStarts", "allSequences"], inplace=True)

        merged = f'%s/pattern_merged_%s/{wildcards.assembly}_%s_words_length_%s.merged.txt' % (out, protocol, protocol, kmer_length)
        df_merged.to_csv(merged, mode="w", index=False, sep=",")
