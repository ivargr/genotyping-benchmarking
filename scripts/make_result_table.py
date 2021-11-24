import logging
import math
logging.basicConfig(level=logging.INFO)
logging.info("Making tables")
from tabulate import tabulate

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
    "usN1000": ["mapI1000", "usN1000"],
    "kage": ["mapI1000", "usN2058"],
    "pangenie": ["pangenieN32"]
}

methods = method_names.split(",")

run_times = {}
memory_usage = {}
accuracy = {}

def is_unix_time_report(file_name):
    with open(file_name) as f:
        content = f.read()
        if "Command being timed:" in content:
            return True

    return False


def get_run_time(job_name, experiment, dataset):
    file_name = "data/" + dataset + "/benchmarks/" + job_name + "_" + experiment + ".tsv"
    if is_unix_time_report(file_name):
        lines = open(file_name).readlines()
        line = [line for line in lines if "Elapsed (wall clock) time (h:mm:ss or m:ss): " in line][0]
        time_string = line.replace("Elapsed (wall clock) time (h:mm:ss or m:ss): ", "").strip().split(":")
        if len(time_string) == 2:
            logging.info(time_string)
            hours = 0
            minutes = int(time_string[0])
            seconds = int(time_string[1].split(".")[0])
        elif len(time_string) == 3:
            hours = int(time_string[0])
            minutes = int(time_string[1])
            seconds = int(time_string[1])
        logging.info("%d hours and %d minutes" % (hours, minutes))
        return hours + minutes/60 + seconds / 3600
    else:
        lines = list(open(file_name).readlines())
        line = lines[1].split()
        run_time = float(line[0]) / (60*60)
    return run_time


def get_memory(job_name, experiment, dataset):
    file_name = "data/" + dataset + "/benchmarks/" + job_name + "_" + experiment + ".tsv"
    if is_unix_time_report(file_name):
        line = [line for line in list(open(file_name).readlines()) if line.startswith("\tMaximum resident set size (kbytes): ")][0]
        memory = float(line.replace("Maximum resident set size (kbytes): ", "").strip()) / 1000000  # gb
    else:
        lines = list(open(file_name).readlines())
        line = lines[1].split()
        memory = round(float(line[2])/1000, 2)  # GBs
    return memory

def get_accuracy(method_name, only_callable_variants=""):
    #file_name = "data/" + dataset + "/happy-" + truth_dataset + "-" + method_name + "_" + experiment + "-only-callable.summary.csv"
    file_name = "data/" + dataset + "/happy-" + truth_dataset + "-" + method_name + "_" + experiment + only_callable_variants + ".summary.csv"
    logging.info("Using file name %s" % file_name)
    lines = list(open(file_name).readlines())
    indels = lines[1].split(",")
    indels_recall = indels[10]
    indels_precision = indels[11]
    indels_f1 = indels[13]
    snps = lines[3].split(",")
    snps_recall = snps[10]
    snps_precision = snps[11]
    snps_f1 = snps[13]
    return [indels_recall, indels_precision, indels_f1, snps_recall, snps_precision, snps_f1]



def format_run_time(hours):
    if hours < 1:
        minutes = math.floor(hours * 60)
        seconds = math.floor(hours * 3600) - minutes*60
        return "%d min" % minutes
    else:
        return "%.1f hours" % hours

def make_table(only_callable_variants=""):
    for method in methods:
        run_times[method] = round(
            sum([get_run_time(job_name, experiment, dataset) for job_name in method_jobs[method]]), 1)
        memory_usage[method] = max([get_memory(job_name, experiment, dataset) for job_name in method_jobs[method]])
        accuracy[method] = get_accuracy(method, only_callable_variants)
        logging.info("Method %s has total run time %d sec and max memory %d bytes" % (
            method, run_times[method], memory_usage[method]))

    table_headers = ["", "Indels recall", "Indels precision", "Indels F1", "SNPs recall", "SNPs precision", "SNPs F1", "Runtime", "Memory usage"]
    table = []

    for method in methods:
        table.append(
            [method.capitalize()] + accuracy[method] + [format_run_time(run_times[method]), str(round(memory_usage[method])) + " GB"]
        )

    print(tabulate(table, table_headers, tablefmt="pretty"))
    print("\n\n\n")
    print(tabulate(table, table_headers, tablefmt="html"))
    print("\n\n\n")
    print(tabulate(table, table_headers, tablefmt="latex"))

make_table()
print("")
make_table("-only-callable")