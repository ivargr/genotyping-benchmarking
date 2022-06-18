import logging
logging.basicConfig(level=logging.INFO)
import sys
n_skipped = 0

for i, line in enumerate(sys.stdin):

    if i % 1000 == 0:
        logging.info("%d lines processed, %d variants with N skipped" % (i, n_skipped))

    if line.startswith("#"):
        print(line, end="")
        continue

    l = line.split()
    if "n" in l[3].lower() or "n" in l[4].lower():
        n_skipped += 1
        continue

    print(line, end="")
