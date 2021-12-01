import logging
logging.basicConfig(level=logging.INFO)
logging.info("Making tables")
from happy_parser import get_accuracy_from_happy_file, get_run_time_from_benchmark_file, get_memory_from_benchmark_file
import sys
import plotly
import plotly.graph_objects as go


method_names = sys.argv[1]
experiment = sys.argv[2]
dataset = sys.argv[3]
truth_dataset = sys.argv[4]
out_file_name = sys.argv[5]


def make_f1_figure(methods, accuracies, file_name):

    colors = ['lightslategray', ] * len(methods)
    x = [m.capitalize() for m in methods]
    y = [accuracies[method]["all"]["f1"] for method in methods]

    fig = go.Figure(data=[go.Bar(
        x=x,
        y=y,
        marker_color=colors  # marker color can be a single color value or an iterable
    )])
    fig.update_layout(title_text='F1 score (2*precision*recall/(precision+recall)')

    plotly.offline.plot(fig, filename=file_name, auto_open=False)
    logging.info("Saved plot to %s" % file_name)


methods = method_names.split(",")

file_names = ["data/" + dataset + "/happy-" + truth_dataset + "-" + method_name + "_" + experiment + ".extended.csv"
              for method_name in methods]
accuracies = {method: get_accuracy_from_happy_file(file_name) for method, file_name in zip(methods, file_names)}

make_f1_figure(methods, accuracies, out_file_name)