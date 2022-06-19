import sys
import functools

# hacky script to sort glimpse chunks so that tabix gets the correct order



def parse_chromosome(chrom):
    if chrom == "X":
        return 23
    elif chrom == "Y":
        return 24
    else:
        return int(chrom)


def sort_function(name1, name2):
    chromosome1 = parse_chromosome(name1.split(".")[-3].split("-")[-1])
    chromosome2 = parse_chromosome(name2.split(".")[-3].split("-")[-1])

    if chromosome1 < chromosome2:
        return -1
    elif chromosome1 == chromosome2:
        id1 = int(name1.split(".")[-2])
        id2 = int(name2.split(".")[-2])
        if id1 < id2:
            return -1

    return 1


files = []
for line in sys.stdin:
    files.append(line.strip())

files.sort(key=functools.cmp_to_key(sort_function))
print("\n".join(files))