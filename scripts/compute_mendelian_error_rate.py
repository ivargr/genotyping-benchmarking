import logging
logging.basicConfig(level=logging.INFO)
import sys
from alignment_free_graph_genotyper.variants import VcfVariants


def read_genotypes_into_dict(file_name, dont_read_genotypes=False):
    logging.info("Reading genotypes from %s" % file_name)
    if dont_read_genotypes:
        genotypes = set()
    else:
        genotypes = {}
    variants = VcfVariants.from_vcf(file_name, make_generator=False, limit_to_chromosome=1, skip_index=True)

    for i, variant in enumerate(variants):
        if i % 100000 == 0:
            logging.info("%d variants processed" % i)
        if dont_read_genotypes:
            genotypes.add(variant.id())
        else:
            genotypes[variant.id()] = variant.get_numeric_genotype()

    return genotypes


all_variants = sys.argv[1]
variants = read_genotypes_into_dict(all_variants, True)
child = read_genotypes_into_dict(sys.argv[2])
mother = read_genotypes_into_dict(sys.argv[3])
father = read_genotypes_into_dict(sys.argv[4])


n_non_homo_ref = 0
n_noncompliant = 0
for variant_id in variants:
    if child[variant_id] != 1:
        n_non_homo_ref += 1
        if child[variant_id] == 2:
            # homo alt
            if mother[variant_id] == 1 or father[variant_id] == 1:
                n_noncompliant += 1
        elif child[variant_id] == 3:
            # heterozygous
            if mother[variant_id] == 1 and father[variant_id] == 1:
                n_noncompliant +=  1
            elif mother[variant_id] == 2 and father[variant_id] == 2:
                n_noncompliant += 1

logging.info("Number of noncompliant variants: %d" % n_noncompliant)
logging.info("Total number of variants that are not homozygous ref: %d" % n_non_homo_ref)