import logging
logging.basicConfig(level=logging.INFO)
import sys
from graph_kmer_index import KmerIndex, ReverseKmerIndex
from obgraph.variant_to_nodes import VariantToNodes
from obgraph.variants import VcfVariants
import pickle


def find_variants_with_nonunique_kmers(variants, variant_to_nodes, population_kmers, variant_kmers):
    nonunique = set()
    for variant_id, variant in enumerate(variants):
        if variant.type != "SNP":
            continue
        ref_node = variant_to_nodes.ref_nodes[variant_id]
        var_node = variant_to_nodes.var_nodes[variant_id]
        kmers = list(variant_kmers.get_node_kmers(ref_node)) + list(variant_kmers.get_node_kmers(var_node))
        kmer_frequencies = [population_kmers.get_frequency(k) for k in kmers]
        if len(kmer_frequencies) > 0 and max(kmer_frequencies) > 1:
            nonunique.add(variant_id)


        if variant_id % 100 == 0:
            logging.info("%d variants checked, %d have nonunique kmers" % (variant_id, len(nonunique)))

    return nonunique



variants = VcfVariants.from_vcf(sys.argv[1])
variant_kmers = ReverseKmerIndex.from_file(sys.argv[2])
population_kmers = KmerIndex.from_file(sys.argv[3])
variant_to_nodes = VariantToNodes.from_file(sys.argv[4])

nonunique = find_variants_with_nonunique_kmers(variants, variant_to_nodes, population_kmers, variant_kmers)
logging.info("Found %d SNPs with nonunique kmers" % len(nonunique))

pickle.dump(nonunique, open(sys.argv[5], "wb"))
logging.info("Wrote to file")
