"""
File is direct copy from Pangenie repository: https://bitbucket.org/jana_ebler/pangenie/raw/d340042099a3af7683729327235a3cd3c8bb8a71/pipelines/run-from-callset/scripts/prepare-vcf.py
"""

import sys
import argparse

parser = argparse.ArgumentParser(prog='prepare-vcf.py', description="cat <vcf-file> | python3 prepare-vcf.py ")
parser.add_argument('--missing', metavar='MISSING', type=float, default=0.0, help="Maximum allowed fraction of missing alleles per position.")
args = parser.parse_args()

total_records = 0
total_alleles = 0
missing_records = 0
missing_alleles = 0
ns_records = 0
ns_alleles = 0
unphased_records = 0
unphased_alleles = 0
written_records = 0
written_alleles = 0

for line in sys.stdin:
	if line.startswith('#'):
		# store header lines
		print(line.strip())
		continue
	total_records += 1
	fields = line.strip().split()
	n_alt_alleles = len(fields[4].split(','))
	total_alleles += n_alt_alleles
	if any([c not in 'CAGTcagt,' for c in fields[4]]):
		ns_records += 1
		ns_alleles += n_alt_alleles
		continue
	n_paths = len(fields[9:]) * 2
	n_missing = 0
	n_total = 0
	for gt in fields[9:]:
		n_total += 2
		if gt == './.':
			n_missing += 2
		elif gt == '.':
			n_missing += 2
		elif '|' in gt:
			alleles = gt.split('|')
			n_missing += alleles.count('.')
		else:
			raise Exception('VCF contains unphased positions.')
	frac_missing = n_missing / n_total
	if frac_missing > args.missing:
		continue
	print(line.strip())
	written_records += 1
	written_alleles += n_alt_alleles

# print statistics
sys.stderr.write('skipped ' + str(missing_records) + ' (' + str(missing_alleles) + ') records (alleles) for which fraction of missing alleles exceeds threshold.')
sys.stderr.write('skipped ' + str(ns_records) + ' (' + str(ns_alleles) + ') records (alleles) for which alternative alleles contained Ns.')
sys.stderr.write('kept ' + str(written_records) + ' (' + str(written_alleles) + ') records (alleles) of ' + str(total_records) + ' (' + str(total_alleles) + ').')
