import logging
logging.basicConfig(level=logging.INFO)
import sys


f = open(sys.argv[1])
regions = sys.argv[2]

n_skipped_lines = 0

if "," in regions:
    logging.info("Assuming regions is a list of chromosomes")
    chromosomes = set(regions.split(","))

    for line in f:
        if line.split()[0] in chromosomes:
            print(line.strip())
        else:
            n_skipped_lines += 1


else:
    logging.info("Assuming regions is a single chromosome region")
    region_chromosome = regions.split(":")[0]
    region_start = int(regions.split(":")[1].split("-")[0])
    region_end = int(regions.split(":")[1].split("-")[1])

    logging.info("Region is chromosome %s from %d to %d" % (region_chromosome, region_start, region_end))

    for line in f:
        l = line.split()
        chromosome = l[0]
        start = int(l[1])
        end = int(l[2])

        if chromosome != region_chromosome:
            n_skipped_lines += 1
            continue

        if start < region_start:
            if end > region_end:
                print(chromosome + "\t" + str(region_start) + "\t" + str(region_end))
            elif end > region_start:
                print(chromosome + "\t" + str(region_start) + "\t" + str(end))
            else:
                n_skipped_lines += 1
            continue
        elif end > region_end:
            # todo fix
            n_skipped_lines += 1
            continue

        print(line.strip())

logging.info("Skipped %d lines in total" % n_skipped_lines)
