import logging
logging.basicConfig(level=logging.INFO)
import sys


max_length = int(sys.argv[1])
n_removed = 0
logging.info("Will filter at length %d" % max_length)

for i, line in enumerate(sys.stdin):
    if i % 10000 == 0:
        logging.info("%d lines processed" % i)

    if line.startswith("#"):
        print(line.strip())
        continue

    l = line.split()
    ref_sequence = l[3]
    alt_sequence = l[4]

    if len(ref_sequence)-1 > max_length or len(alt_sequence)-1 > max_length:
        n_removed += 1
        continue

    print(line.strip())

logging.info("%d indels removed" % n_removed)