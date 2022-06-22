import sys
from graph_kmer_index.index_bundle import IndexBundle
bundle = IndexBundle.from_file(snakemake.input[0])
bundle.indexes["KmerIndex"] = None
bundle.to_file(snakemake.output[0], compress=False)
