import sys
import random

file = sys.argv[1]
seed = int(sys.argv[2])
random.seed(seed)

with open(file) as f:
    lines = list(f.readlines())
    random.shuffle(lines)

    print('\n'.join((line.strip() for line in lines)))
