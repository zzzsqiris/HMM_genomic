import HMM_utils
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("pred_file")
parser.add_argument("gff_file")
args = parser.parse_args()

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

    for p, t in zip(predicted_path, true_path):
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

    sn = tp / (tp + fn)
    sp = tp / (tp + fp)
    acc = (tp + tn) / (tp + tn + fp + fn)

    return {
        "Mode": "Multi" if is_multi_state else "2-state",
        "Counts": {"TP": tp, "TN": tn, "FP": fp, "FN": fn},
        "Sn": round(sn, 4),
        "Sp": round(sp, 4),
        "Acc": round(acc, 4)
    }

with open(args.pred_file, 'r') as f:
        content = f.read().strip()
        predicted_path = content.split('\t')

results = evaluate_performance(predicted_path, args.gff_file)
print(results)