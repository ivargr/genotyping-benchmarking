import sys
import logging

rate = int(sys.argv[1])


i = 0
skip = True
for line_number, line in enumerate(sys.stdin):

    if line_number % 4 == 0:
        if i % rate == 0:
            skip = False         
        else:
            skip = True

        i += 1

    if not skip:
        print(line.strip())
