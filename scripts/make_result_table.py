import logging
import sys


method_names = sys.argv[1]
experiment = sys.argv[2]
dataset = sys.argv[3]
truth_dataset = sys.argv[4]

# Mapping from method name to jobs for that method
method_jobs  = {
    "bayestyper": ["bayestyper_kmc", "bayestyper_bloomfilter", "bayestyper"],
    "gatk": ["bwamem", "gatk"],
    "graphtyper": ["bwamem", "graphtyper"],
    "malva": ["malva_kmc", "malva"],
    "us": ["mapI1000", "usN2058"],
    "pangenie": ["pangenieN32"]
}

methods = method_names.split(",")

run_times = {}
memory_usage = {}
accuracy = {}

def get_run_time(job_name, experiment, dataset):
    file_name = "data/" + dataset + "/benchmarks/" + job_name + "_" + experiment + ".tsv"
    lines = list(open(file_name).readlines())
    line = lines[1].split()
    run_time = float(line[0]) / (60*60)
    return run_time


def get_memory(job_name, experiment, dataset):
    file_name = "data/" + dataset + "/benchmarks/" + job_name + "_" + experiment + ".tsv"
    lines = list(open(file_name).readlines())
    line = lines[1].split()
    memory = round(float(line[2])/1000, 2)  # GBs
    return memory

def get_accuracy(method_name):
    file_name = "data/" + dataset + "/happy-" + truth_dataset + "-" + method_name + "_" + experiment + ".summary.csv"
    lines = list(open(file_name).readlines())
    indels = lines[1].split(",")
    indels_recall = indels[10]
    indels_precision = indels[11]
    snps = lines[3].split(",")
    snps_recall = snps[10]
    snps_precision = snps[11]
    return [indels_recall, indels_precision, snps_recall, snps_precision]

for method in methods:
    run_times[method] = round(sum([get_run_time(job_name, experiment, dataset) for job_name in method_jobs[method]]), 1)
    memory_usage[method] = max([get_memory(job_name, experiment, dataset) for job_name in method_jobs[method]])
    accuracy[method] = get_accuracy(method)
    logging.info("Method %s has total run time %d sec and max memory %d bytes" % (method, run_times[method], memory_usage[method]))


table_headers = ["", "Indels recall", "Indels precision", "SNPs recall", "SNPs precision", "Runtime", "Memory usage"]
table = []

for method in methods:
    table.append(
        [method.capitalize()] + accuracy[method] + [str(run_times[method]) + " h", str(memory_usage[method]) + " GB"]
    )

from tabulate import tabulate

print(tabulate(table, table_headers, tablefmt="pretty"))
print("\n\n\n")
print(tabulate(table, table_headers, tablefmt="html"))
print("\n\n\n")
print(tabulate(table, table_headers, tablefmt="latex"))
