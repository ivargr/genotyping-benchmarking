import sys

for line in sys.stdin:
    if not line.startswith(">"):
        print(line.strip())
        continue

    print(line.split(":")[0].strip())
