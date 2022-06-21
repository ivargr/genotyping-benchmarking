import logging
logging.basicConfig(level=logging.INFO)
import sys
from graph_kmer_index import KmerIndex, ReverseKmerIndex
from obgraph.variant_to_nodes import VariantToNodes
from obgraph.variants import VcfVariants
import pickle
from tabulate import tabulate


def get_variants_with_matching_genotype(nonuniqe, all_variants, truth, query_variants):
    n_total_unique = 0
    n_total_nonunique = 0
    n_correct = 0
    n_correct_unique = 0
    n_correct_nonunique = 0
    for i, variant in enumerate(all_variants):
        if truth.has_variant(variant):
            if query_variants.has_variant_genotype(truth.get(variant)):
                n_correct += 1

                if i in nonuniqe:
                    n_correct_nonunique += 1
                else:
                    n_correct_unique += 1

            if i in nonuniqe:
                n_total_nonunique += 1
            else:
                n_total_unique += 1

    return n_total_unique, n_total_nonunique, n_correct_unique, n_correct_nonunique




variants = VcfVariants.from_vcf(sys.argv[1])
naive_kage = VcfVariants.from_vcf(sys.argv[2])
kage_modelled_counts = VcfVariants.from_vcf(sys.argv[3])
kage = VcfVariants.from_vcf(sys.argv[4])
truth = VcfVariants.from_vcf(sys.argv[5])
nonunique = pickle.load(open(sys.argv[6], "rb"))


logging.info("Found %d SNPs with nonunique kmers" % len(nonunique))

stats = [
    (name, get_variants_with_matching_genotype(nonunique, variants, truth, k))
    for name, k in
         [
            ("NAIVE KAGE", naive_kage),
            ("KAGE with model of kmer counts", kage_modelled_counts),
            ("KAGE (full)", kage)
        ]
]

table = []
table_headers = ["Method", "Accuracy on SNPs with non-unique kmers", "Accuracy on SNPs with unique kmers"]

for name, s in stats:
    logging.info("--- %s ---" % name)
    accuracy_unique = 100 * s[2] / s[0]
    accuracy_nonunique = 100 * s[3] / s[1]
    logging.info("Accuracy on nonunique: %.3f" % (accuracy_nonunique))
    logging.info("Accuracy on unique: %.3f" % (accuracy_unique))
    table.append([name, str(round(accuracy_nonunique)) + " % ", str(round(accuracy_unique)) + " % "])


    print(tabulate(table, table_headers, tablefmt="pretty"))
    print("\n\n\n")
    print(tabulate(table, table_headers, tablefmt="html"))
    print("\n\n\n")
    print(tabulate(table, table_headers, tablefmt="latex"))




