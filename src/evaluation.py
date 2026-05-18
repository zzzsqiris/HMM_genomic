import HMM_utils
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("pred_file")
parser.add_argument("gff_file")
args = parser.parse_args()

def read_pred_gff(pred_file):
    pred_info = []

    with open(pred_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line == "":
                continue

            parts = line.split('\t')
            ftype = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            pred_info.append((ftype, start, end))

    seq_length = 0
    for ftype, start, end in pred_info:
        if end > seq_length:
            seq_length = end

    predicted_path = ['Intron'] * seq_length

    for ftype, start, end in pred_info:
        for i in range(start - 1, end):
            predicted_path[i] = ftype

    return predicted_path

def evaluate_performance(predicted_path, gff3_file_path):
    seq_length = len(predicted_path)

    # determine simple or multiple states
    states_set = set(predicted_path)
    if len(states_set) > 2:
        is_multi_state = True
    else:
        is_multi_state = False

    # build true path from gff3
    true_path = ['Intron'] * seq_length
    gff_into = HMM_utils.read_gff(gff3_file_path)

    for ftype, start, end in gff_into:
        target_types = ['exon', 'five_prime_UTR', 'three_prime_UTR', 'CDS']
        if ftype in target_types:
            if is_multi_state:
                label = ftype 
            else:
                label = 'Exon'
            for i in range(start - 1, end):
                if i < seq_length:
                    true_path[i] = label


# compare True path and Predicted path

#                  True path
#                E         I
#   Pred   E     TP       FP  
#   path   I     FN       TN    

    tp, tn, fp, fn = 0, 0, 0, 0

    for i in range(seq_length):
        p = predicted_path[i]
        t = true_path[i]
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

    return {
        "Mode": "Multi" if is_multi_state else "2-state",
        "Counts": {"TP": tp, "TN": tn, "FP": fp, "FN": fn},
        "Sn": round(sn, 4),
        "Sp": round(sp, 4),
        "Acc": round(acc, 4)
    }

predicted_path = read_pred_gff(args.pred_file)

results = evaluate_performance(predicted_path, args.gff_file)
# print(results)

out_file = args.pred_file.replace(".pred.gff3", ".eval.txt")

with open(out_file, "w") as f:
    counts = results["Counts"]
    f.write("Sn\tSp\tAcc\tTP\tTN\tFP\tFN\n")
    f.write(str(results["Sn"]) + "\t")
    f.write(str(results["Sp"]) + "\t")
    f.write(str(results["Acc"]) + "\t")
    f.write(str(counts["TP"]) + "\t")
    f.write(str(counts["TN"]) + "\t")
    f.write(str(counts["FP"]) + "\t")
    f.write(str(counts["FN"]) + "\n")
