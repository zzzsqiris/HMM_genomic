import numpy
import math
import HMM_utils
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("json_prob_file")
parser.add_argument("fa_file")
parser.add_argument("--out_dir", default="build")
args = parser.parse_args()

prob_file = args.json_prob_file
input_path = args.fa_file
out_dir = args.out_dir

states, TP, EP = HMM_utils.read_in_prob(prob_file)
num_states = len(states)

nt_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
init_log = [math.log(1 / num_states)] * num_states

def predict_one_file(fa_file):
    fasta_dict = HMM_utils.read_fasta(fa_file)
    raw_seq = list(fasta_dict.values())[0].upper()
    seq = [nt_map[nt] for nt in raw_seq if nt in nt_map]

    dpm = numpy.zeros((num_states, len(seq) + 1))
    trace = numpy.full((num_states, len(seq) + 1), -1)

    for s in range(num_states):
        dpm[s][0] = init_log[s]

    # loop through seq
    for i in range(1, len(seq)+1): 
        obs = seq[i-1]
        # loop through states
        for cur in range(num_states):
            # init as first state
            max_log_p = dpm[0][i-1] + TP[0][cur] + EP[cur][obs]
            best_prev_node = 0
            
            # loop through other states
            for prev in range(1, num_states):
                current_p = dpm[prev][i-1] + TP[prev][cur] + EP[cur][obs]
                if current_p > max_log_p:
                    max_log_p = current_p
                    best_prev_node = prev
            
            dpm[cur][i] = max_log_p
            trace[cur][i] = best_prev_node

    # find last state, init as first state
    max_val = dpm[0][len(seq)]
    current = 0

    for s in range(1, num_states):
        if dpm[s][len(seq)] > max_val:
            max_val = dpm[s][len(seq)]
            current = s

    path = [current]

    # trace back
    for i in range(len(seq), 0, -1):
        prev_state = trace[current][i]
        if prev_state != -1:
            path.append(prev_state)
            current = prev_state

    path.reverse()
    path_names = []
    for p in path[1:]:
        name = states[p]
        name = HMM_utils.state_for_output(name)
        path_names.append(name)

    seq_id = list(fasta_dict.keys())[0]

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    out_file = os.path.join(out_dir, seq_id + ".pred.gff3")

    with open(out_file, "w") as f:
        start = 1
        current_state = path_names[0]

        for i in range(1, len(path_names)):
            if path_names[i] != current_state:
                line = [seq_id, 'HMM_pred', current_state, str(start), str(i), '.', '.', '.', '.']
                f.write('\t'.join(line) + '\n')
                start = i + 1
                current_state = path_names[i]

        line = [seq_id, 'HMM_pred', current_state, str(start), str(len(path_names)), '.', '.', '.', '.']
        f.write('\t'.join(line) + '\n')

if os.path.isdir(input_path):
    for filename in os.listdir(input_path):
        if filename.endswith(".fa"):
            fa_file = os.path.join(input_path, filename)
            predict_one_file(fa_file)
else:
    predict_one_file(input_path)
