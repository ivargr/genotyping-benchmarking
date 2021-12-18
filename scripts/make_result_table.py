import logging
import math

import plotly

logging.basicConfig(level=logging.INFO)
logging.info("Making tables")
from tabulate import tabulate
from happy_parser import get_accuracy_from_happy_file, get_run_time_from_benchmark_file, get_memory_from_benchmark_file
import sys
import plotly.express as px


method_names = sys.argv[1]
experiment = sys.argv[2]
dataset = sys.argv[3]
truth_dataset = sys.argv[4]
out_base_name = sys.argv[5]


name_mappings = {"bayestyper": "Bayestyper",
                 "gatk": "GATK",
                 "us": "KAGE",
                 "pangenie": "PanGenie",
                 "malva": "Malva",
                 "graphtyper": "Graphtyper"}

# Mapping from method name to jobs for that method
method_jobs  = {
    "bayestyper": ["bayestyper_kmc", "bayestyper_bloomfilter", "bayestyper"],
    "gatk": ["bwamem", "gatk"],
    "graphtyper": ["bwamem", "graphtyper"],
    "malva": ["malva_kmc", "malva"],
    "us": ["mapI1000", "usN2548all"],
    "usN1000": ["mapI1000", "usN1000"],
    "KAGE": ["mapI1000", "usN2548all"],
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



def make_plots():

    accuracies = []
    memories = []
    run_times = []
    for method in methods:
        file_name = "data/" + dataset + "/happy-" + truth_dataset + "-" + method + "_" + experiment + ".extended.csv"
        r = get_accuracy_from_happy_file(file_name)
        accuracies.append(r["all"]["f1"])
        memories.append(max([get_memory_from_benchmark_file(job_name, experiment, dataset) for job_name in method_jobs[method]]))
        run_times.append(round(
            sum([get_run_time_from_benchmark_file(job_name, experiment, dataset) for job_name in method_jobs[method]]), 2))

    method_names = [name_mappings[name] for name in methods]
    for plot_type, yaxis_name, data in zip(["f1", "memory", "runtime"],
                                           ["F1 score", "Max memory (GB)", "Runtime (hours)"],
                                           [accuracies, memories, run_times]):

        fig = px.bar(x=method_names, y=data)
        fig.update_layout(
            yaxis=dict(
                tickfont=dict(size=20),
                showgrid=False
            ),
            xaxis=dict(showgrid=False),
            font=dict(size=20),
            xaxis_title="Method",
            yaxis_title=yaxis_name,
        )

        fig.update_xaxes(showline=True, linewidth=1, linecolor='#666666')
        fig.update_yaxes(showline=True, linewidth=2, linecolor='#666666')

        fig.update_traces(marker_color='lightslategray')
        if plot_type == "f1":
            fig.update_yaxes(range=[0.8, 1])
        plotly.offline.plot(fig, filename=plot_type + out_base_name, auto_open=False)



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
        table.append(
            [name_mappings[method]] + accuracy[method] + [format_run_time(run_times[method]), str(round(memory_usage[method])) + " GB"]
        )

    print(tabulate(table, table_headers, tablefmt="pretty"))
    print("\n\n\n")
    print(tabulate(table, table_headers, tablefmt="html"))
    print("\n\n\n")
    print(tabulate(table, table_headers, tablefmt="latex"))


make_table()
#print("<h3>Only including callable variants</h3>")
#make_table("-only-callable")
#print("<h3>Long indels not included</h3>")
#make_table("-short-indels")

make_plots()