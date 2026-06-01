import os
import argparse
import HMM_utils

parser = argparse.ArgumentParser()
parser.add_argument("data_dir")
parser.add_argument("--tolerance", type=int, default=5)
parser.add_argument("--pred_dir", default="build")
args = parser.parse_args()

data_dir = args.data_dir
tolerance = args.tolerance
pred_dir = args.pred_dir

def find_switches(path):
    switches = []

    for i in range(1, len(path)):
        if path[i] == 'Skip' or path[i - 1] == 'Skip':
            continue

        if path[i] != path[i - 1]:
            switches.append(i + 1)

    return switches

def count_matched(true_switches, pred_switches, tolerance):
    matched = 0
    used_pred = []

    for true_pos in true_switches:
        for pred_pos in pred_switches:
            if pred_pos in used_pred:
                continue

            if abs(true_pos - pred_pos) <= tolerance:
                matched += 1
                used_pred.append(pred_pos)
                break

    return matched

sensitivity_sum = 0
specificity_sum = 0
average_sum = 0
num_files = 0
missing_files = 0

for filename in os.listdir(data_dir):
    if filename.endswith(".fa"):
        base_name = filename.replace(".fa", "")

        pred_file = os.path.join(pred_dir, base_name + ".pred.gff3")
        gff_file = os.path.join(data_dir, base_name + ".gff3")

        if not os.path.exists(gff_file):
            continue

        if not os.path.exists(pred_file):
            missing_files += 1
            continue

        predicted_path = HMM_utils.read_pred_gff(pred_file)
        seq_length = len(predicted_path)
        true_path = HMM_utils.read_true_path(gff_file, seq_length, False)

        pred_switches = find_switches(predicted_path)
        true_switches = find_switches(true_path)

        matched = count_matched(true_switches, pred_switches, tolerance)

        if len(true_switches) == 0:
            switch_sn = 0
        else:
            switch_sn = matched / len(true_switches)

        if len(pred_switches) == 0:
            switch_sp = 0
        else:
            switch_sp = matched / len(pred_switches)

        switch_ave = (switch_sn + switch_sp) / 2

        sensitivity_sum += switch_sn
        specificity_sum += switch_sp
        average_sum += switch_ave
        num_files += 1

if num_files == 0:
    print("No files evaluated")
else:
    print("Number of files:", num_files)
    print("Missing prediction files:", missing_files)
    print("Tolerance:", tolerance)
    print("Average switch evaluation")
    print("Switch Sn:", round(sensitivity_sum / num_files, 4))
    print("Switch Sp:", round(specificity_sum / num_files, 4))
    print("Switch Ave:", round(average_sum / num_files, 4))
