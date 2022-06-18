import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')
import argparse
import sys
import glob
import plotly.express as px
import plotly
import plotly.graph_objects as go
import numpy as np
from happy_parser import get_accuracy_from_happy_file



def get_snp_recall_and_precision_from_happy_csv(file_name):
    f = open(file_name)
    lines = f.readlines()
    line = [line for line in lines if line.startswith("SNP,ALL")]
    assert len(line) == 1
    line = line[0]
    print(line)
    line = line.split(",")
    recall = float(line[10])
    precision = float(line[11])
    return recall, precision


def plot_result_files(args):
    logging.info("Plotting")

    labels = []

    color_mappings = {"pangenie": "#8F3333",
                      "us": "#456F9E",
                      "kage": "#456F9E",
                      "kageNoHelperModel": "#456F9E",
                      "malva": "#142321",
                      "naivekage": "#A8D4D4",
                      "model": "#456F9E",
                      "graphtyper": "purple",
                      "bayestyper": "green",
                      "gatk": "gray"}

    name_mappings = {"pangenie": "Pangenie", "us": "KAGE", "naivekage": "Naive KAGE", "kageNoHelperModel": "KAGE", "kage": "KAGE",
                     "model": "Modelled kmer counts", "malva": "Malva", "graphtyper": "Graphtyper",
                     "bayestyper": "Bayestyper", "gatk": "GATK"}


    fig = go.Figure(
        layout=go.Layout(
            xaxis=dict(showgrid=True, zeroline=False),
                #plot_bgcolor='rgba(0, 0, 0, 0)',
                #paper_bgcolor='rgba(235,242,247,1)'
    )
    )

    fig.update_xaxes(showline=True, linewidth=2, linecolor='#666666')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='#666666')

    if args.type == "f1":
        fig.update_layout(xaxis=dict(tickmode="linear"))

    # fig.update_layout(title_text="Reads (%s)" % self.type)
    fig.update_layout(
        xaxis=dict(
            #showexponent='all',
            #exponentformat='e',
            nticks=5,
            tickfont=dict(
                size=20
            ),
            showgrid=False
        )
    )

    if args.type == "f1":
        xaxis = "Number of individuals"
        yaxis = "F1 score"
    else:
        xaxis = "Recall"
        yaxis = "Precision"

    fig.update_layout(
        yaxis=dict(
            tickfont=dict(
                size=20
            ),
            showgrid=False
        ),
        font=dict(size=20),
        xaxis_title=xaxis,
        yaxis_title=yaxis
    )

    if args.type == "f1":
        fig.update_xaxes(type='category')

    fig.update_layout(showlegend=True)
    #fig.update_layout(plot_bgcolor='rgba(230,230,230,255)')
    #fig.update_layout(xaxis_type="log")

    file_names = args.results_file_prefixes.split(",")
    method_names = args.method_names.split(",")
    legends_shown = set()

    max_y = 0
    min_y = 1
    #for prefixes in args.results_file_prefixes.split(","):
    for file_name, method_name in zip(file_names, method_names):
        x = []
        y = []
        colors = []
        ticker_texts = []
        #files = glob.glob(prefixes + "*.txt")

        #logging.info("Will plot %s" % files)
        data = {}
        #for file in files:
            #data = f.readline().split()
            #recall = float(data[0])
            #precision = float(data[1])

        #recall, precision = get_snp_recall_and_precision_from_happy_csv(file_name)
        stats = get_accuracy_from_happy_file(file_name)
        logging.info("Stats for %s: %s" % (method_name, stats))
        recall = stats["all"]["recall"]
        precision = stats["all"]["precision"]
        f1 = stats["all"]["f1"]

        #n_individuals = int(file.split("_")[2].split(".")[0])
        ticker_text = name_mappings[method_name]
        ticker_text = ""

        #if method_name != "nomodel" and method_name != "malva" and method_name != "model":
        #    ticker_text = str(n_individuals)
        n_individuals = 2548
        if "N" in file_name:
            logging.info("File name: %s" % file_name)
            n_individuals = int(file_name.split("-")[2].split("_")[0].split("N")[-1].replace("all", ""))

        ticker_texts.append(ticker_text)

        if args.type == "f1":
            x.append(n_individuals)
            y.append(f1)

        else:
            x.append(recall)
            y.append(precision)


        #colors.append(color_mappings[file.split("_")[1]])
        chosen_color = "black"
        for possible_color_name, color in color_mappings.items():
            if possible_color_name in method_name:
                chosen_color = color
                break

        logging.info("Method name: " + method_name)

        if "nomodel" in method_name:
            ticker_text = "Naive method"
            chosen_color = "orange"

        colors.append(chosen_color)

        name = name_mappings[method_name]

        showlegend = True
        if name in legends_shown:
            showlegend = False

        legends_shown.add(name)


        fig.add_trace(go.Scatter(x=x, y=y,
                                 mode='markers',
                                 textposition="bottom left",
                                 name=name,
                                 showlegend=showlegend,
                                 textfont=dict(
                                     size=15
                                 ),
                                 marker=dict(
                                     size=20,
                                     color=color_mappings[method_name]
                                 )
                        ),
                      )


        min_y = min([min_y, min(y)])
        max_y = max([max_y, max(y)])

    if args.type == "f1":
        logging.info("Max y value is %.4f" % max_y)
        fig.update_yaxes(range=[min_y - 0.01, max_y+0.01])
    else:
        fig.update_yaxes(range=[min_y-0.005, 1.0])
    save_to_file = args.out_file_name
    if save_to_file is not None:
        plotly.offline.plot(fig, filename=save_to_file, auto_open=False)
        logging.info("Saved plot to %s" % save_to_file)
    else:
        fig.show()



def run_argument_parser(args):
    parser = argparse.ArgumentParser(
        description='Alignment free graph genotyper',
        prog='alignment_free_graph_genotyper',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=50, width=100))

    subparsers = parser.add_subparsers()
    subparser = subparsers.add_parser("plot_results_files")
    subparser.add_argument("-f", "--results_file_prefixes", required=True, help="Comma-separated list of result file names")
    subparser.add_argument("-n", "--method-names", required=True, help="Comma-separated list of method names")
    subparser.add_argument("-o", "--out-file-name")
    subparser.add_argument("-t", "--type", default="scatter", help="scatter or f1")
    subparser.set_defaults(func=plot_result_files)


    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(args)
    args.func(args)

run_argument_parser(sys.argv[1:])
