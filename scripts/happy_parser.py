import logging


def get_accuracy_from_happy_file(happy_csv_file_name):
    f = open(happy_csv_file_name)
    lines = f.readlines()
    lines = [line.split(",") for line in lines if line.startswith("SNP,*,*,ALL") or line.startswith("INDEL,*,*,ALL")]
    assert len(lines) == 2

    query_tp_sum = 0
    query_fp_sum = 0
    truth_tp_sum = 0
    truth_fn_sum = 0
    results = {}
    for line in lines:
        type = line[0].lower()
        truth_tp = float(line[23])
        truth_fn = float(line[30])
        query_tp = float(line[44])
        query_fp = float(line[51])

        truth_tp_sum += truth_tp
        truth_fn_sum += truth_fn
        query_tp_sum += query_tp
        query_fp_sum += query_fp

        local_precision = query_tp / (query_tp + query_fp)
        local_recall = truth_tp / (truth_tp + truth_fn)
        local_f1 = 2*local_precision*local_recall / (local_recall+local_precision)
        results[type] = {"precision": local_precision, "recall": local_recall, "f1": local_f1}

    global_precision = query_tp_sum / (query_tp_sum + query_fp_sum)
    global_recall = truth_tp_sum / (truth_tp_sum + truth_fn_sum)
    global_f1 = 2*global_recall*global_precision/(global_recall+global_precision)

    results["all"] = {"precision": global_precision, "recall": global_recall, "f1": global_f1}
    return results


def is_unix_time_report(file_name):
    with open(file_name) as f:
        content = f.read()
        if "Command being timed:" in content:
            return True

    return False


def get_run_time_from_benchmark_file(job_name, experiment, dataset):
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
        logging.info("%d hours and %d minutes and %d seconds (%s)" % (hours, minutes, seconds, time_string))
        time_in_hours = hours + minutes/60 + seconds / 3600
        logging.info("Time in hours: %.4f" % time_in_hours)
        return time_in_hours
    else:
        lines = list(open(file_name).readlines())
        line = lines[1].split()
        run_time = float(line[0]) / (60*60)
    return run_time


def get_memory_from_benchmark_file(job_name, experiment, dataset):
    file_name = "data/" + dataset + "/benchmarks/" + job_name + "_" + experiment + ".tsv"
    if is_unix_time_report(file_name):
        line = [line for line in list(open(file_name).readlines()) if line.startswith("\tMaximum resident set size (kbytes): ")][0]
        memory = float(line.replace("Maximum resident set size (kbytes): ", "").strip()) / 1000000  # gb
    else:
        lines = list(open(file_name).readlines())
        line = lines[1].split()
        memory = round(float(line[2])/1000, 2)  # GBs
    return memory

