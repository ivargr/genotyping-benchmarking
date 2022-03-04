import logging
logging.basicConfig(level=logging.INFO)
import sys


n_skipped = 0
n_kept = 0

for line in sys.stdin:
    if line.startswith("#"):
        print(line.strip())
        continue

    if ">" in line:
        n_skipped += 1
        continue

    print(line.strip())
    n_kept += 1


logging.info("%d skipped" %  n_skipped)
logging.info("%d kept" %  n_kept)
