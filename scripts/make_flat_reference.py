import logging
logging.basicConfig(level=logging.INFO)
import sys


f = open(sys.argv[1])

logging.info("Reading lines")
lines = (line.strip() for line in f if not line.startswith(">"))

print(">ref")
for line in lines:
    print(line, end="")