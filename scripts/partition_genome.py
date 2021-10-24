import sys                                                                                                                                                                                                                                                                                
from collections import OrderedDict
import numpy as np
import math

chunk_size = int(sys.argv[1])
fai = sys.argv[2]

wrapper_text = None
if len(sys.argv) > 3:
    wrapper_text = sys.argv[3]


chromosomes = OrderedDict()

with open(fai) as f:
    for line in f:
        l = line.split()
        chromosome = l[0]
        size = int(l[1])

        chromosomes[chromosome] = size


for chromosome, size in chromosomes.items():
    n_boundaries = max(2, math.ceil(size/chunk_size))
    boundaries = [int(n) for n in np.linspace(1, size, n_boundaries)]
    intervals = [(start, end) for start, end in zip(boundaries[0:-1], boundaries[1:])]
    for interval in intervals:
        interval = "%s:%d-%d" % (chromosome, interval[0], interval[1])
        if wrapper_text is not None:
            print(wrapper_text.replace("---", interval))
        else:
            print(interval)