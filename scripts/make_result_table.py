import logging
import math
logging.basicConfig(level=logging.INFO)
logging.info("Making tables")
from tabulate import tabulate
from happy_parser import get_accuracy_from_happy_file, get_run_time_from_benchmark_file, get_memory_from_benchmark_file
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
    "KAGE": ["mapI1000", "usN2058"],
    "pangenie": ["pangenieN32"]
}

methods = method_names.split(",")


def get_accuracy(method_name, only_callable_variants=""):
    file_name = "data/" + dataset + "/happy-" + truth_dataset + "-" + method_name + "_" + experiment + only_callable_variants + ".extended.csv"
    logging.info("Using file name %s" % file_name)
    r = get_accuracy_from_happy_file(file_name)
    results = [r["indel"]["recall"], r["indel"]["precision"], r["indel"]["f1"], r["snp"]["recall"], r["snp"]["precision"], r["snp"]["f1"]]
    results = ["%.3f" % float(r) for r in results]
    logging.info("Accuracy: %s" % results)
    return results




def format_run_time(hours):
    if hours < 1:
        minutes = math.floor(hours * 60)
        seconds = math.floor(hours * 3600) - minutes*60
        return "%d min" % minutes
    else:
        return "%.1f hours" % hours


def make_table(only_callable_variants=""):
    run_times = {}
    memory_usage = {}
    accuracy = {}
    for method in methods:
        run_times[method] = round(
            sum([get_run_time_from_benchmark_file(job_name, experiment, dataset) for job_name in method_jobs[method]]), 2)
        memory_usage[method] = max([get_memory_from_benchmark_file(job_name, experiment, dataset) for job_name in method_jobs[method]])
        accuracy[method] = get_accuracy(method, only_callable_variants)
        logging.info("Method %s has total run time %d sec and max memory %d bytes and accuracy %s" % (
            method, run_times[method], memory_usage[method], accuracy[method]))

    table_headers = ["", "Indels recall", "Indels precision", "Indels F1", "SNPs recall", "SNPs precision", "SNPs F1", "Runtime", "Memory usage"]
    table = []

    for method in methods:
        method_name = method
        if method_name == "us":
            method_name = "KAGE"
        table.append(
            [method_name.capitalize()] + accuracy[method] + [format_run_time(run_times[method]), str(round(memory_usage[method])) + " GB"]
        )

    print(tabulate(table, table_headers, tablefmt="pretty"))
    print("\n\n\n")
    print(tabulate(table, table_headers, tablefmt="html"))
    print("\n\n\n")
    print(tabulate(table, table_headers, tablefmt="latex"))


make_table()
print("<h3>Only including callable variants</h3>")
make_table("-only-callable")
print("<h3>Long indels not included</h3>")
make_table("-short-indels")
