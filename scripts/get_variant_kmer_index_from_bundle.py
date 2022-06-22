from graph_kmer_index.index_bundle import IndexBundle
bundle = IndexBundle.from_file(snakemake.input[0])
kmer_index = bundle.indexes["KmerIndex"]
kmer_index.to_file(snakemake.output[0])
