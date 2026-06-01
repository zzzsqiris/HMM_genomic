import os
import argparse
import HMM_utils

parser = argparse.ArgumentParser()
parser.add_argument("data_dir")
parser.add_argument("--pred_dir", default="build")
args = parser.parse_args()

data_dir = args.data_dir
pred_dir = args.pred_dir

def evaluate_performance(predicted_path, gff_file):
    seq_length = len(predicted_path)

    states_set = set(predicted_path)
    if len(states_set) > 2:
        is_multi_state = True
    else:
        is_multi_state = False

    true_path = HMM_utils.read_true_path(gff_file, seq_length, is_multi_state)

    tp, tn, fp, fn = 0, 0, 0, 0

    for i in range(seq_length):
        p = predicted_path[i]
        t = true_path[i]

        if t == 'Skip':
            continue

        if p == t:
            if p != 'Intron':
                tp += 1
            else:
                tn += 1
        else:
            if p != 'Intron' and t == 'Intron':
                fp += 1
            elif p == 'Intron' and t != 'Intron':
                fn += 1
            else:
                fp += 1

    if tp + fn == 0:
        sn = 0
    else:
        sn = tp / (tp + fn)

    if tp + fp == 0:
        sp = 0
    else:
        sp = tp / (tp + fp)

    if tp + tn + fp + fn == 0:
        acc = 0
    else:
        acc = (tp + tn) / (tp + tn + fp + fn)

    return sn, sp, acc

sensitivity_sum = 0
specificity_sum = 0
acc_sum = 0
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
        sn, sp, acc = evaluate_performance(predicted_path, gff_file)

        sensitivity_sum += sn
        specificity_sum += sp
        acc_sum += acc
        num_files += 1

if num_files == 0:
    print("No files evaluated")
else:
    ave_sn = sensitivity_sum / num_files
    ave_sp = specificity_sum / num_files
    ave_acc = acc_sum / num_files

    print("Number of files:", num_files)
    print("Missing prediction files:", missing_files)
    print("Average")
    print("Sn:", round(ave_sn, 4))
    print("Sp:", round(ave_sp, 4))
    print("Acc:", round(ave_acc, 4))
