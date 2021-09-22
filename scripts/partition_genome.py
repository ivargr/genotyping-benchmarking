import sys                                                                                                                                                                                                                                                                                
from collections import OrderedDict
import numpy as np

chunk_size = int(sys.argv[1])
fai = sys.argv[2]

chromosomes = OrderedDict()

with open(fai) as f:
    for line in f:
        l = line.split()
        chromosome = l[0]
        size = int(l[1])

        chromosomes[chromosome] = size


for chromosome, size in chromosomes.items():
    boundaries = [int(n) for n in np.linspace(1, size, size//chunk_size)]
    intervals = [(start, end) for start, end in zip(boundaries[0:-1], boundaries[1:])]
    for interval in intervals:
        print("%s:%d-%d" % (chromosome, interval[0], interval[1]))
