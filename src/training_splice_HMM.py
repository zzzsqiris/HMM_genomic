import HMM_utils
import os
import json
import math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("data_dir")
args = parser.parse_args()

data_dir = args.data_dir

states = ["Exon", "Donor_G", "Donor_T", "Intron", "Acceptor_A", "Acceptor_G"]
nt_order = ['A', 'C', 'G', 'T']

nt_counts = {}
state_change_counts = {}

for state in states:
    nt_counts[state] = {'A': 1, 'C': 1, 'G': 1, 'T': 1}
    state_change_counts[state] = {}
    for state2 in states:
        state_change_counts[state][state2] = 0

next_states = {
    "Exon": ["Exon", "Donor_G"],
    "Donor_G": ["Donor_T"],
    "Donor_T": ["Intron"],
    "Intron": ["Intron", "Acceptor_A"],
    "Acceptor_A": ["Acceptor_G"],
    "Acceptor_G": ["Exon"]
}

for state in next_states:
    for state2 in next_states[state]:
        state_change_counts[state][state2] = 1

for filename in os.listdir(data_dir):
    if filename.endswith(".fa"):
        fa_path = os.path.join(data_dir, filename)
        gff_path = os.path.join(data_dir, filename[:-3] + ".gff3")

        if not os.path.exists(gff_path):
            continue

        fasta_dict = HMM_utils.read_fasta(fa_path)
        gene_name = filename[:-3]
        dna_seq = fasta_dict[gene_name].upper()

        true_path = HMM_utils.read_true_path(gff_path, len(dna_seq), False)
        features = HMM_utils.read_gff(gff_path)

        for ftype, start, end in features:
            if ftype == "intron":
                if end - start + 1 >= 4:
                    donor_g = start - 1
                    donor_t = start
                    acceptor_a = end - 2
                    acceptor_g = end - 1

                    if donor_g < len(true_path) and true_path[donor_g] == "Intron":
                        true_path[donor_g] = "Donor_G"
                    if donor_t < len(true_path) and true_path[donor_t] == "Intron":
                        true_path[donor_t] = "Donor_T"
                    if acceptor_a < len(true_path) and true_path[acceptor_a] == "Intron":
                        true_path[acceptor_a] = "Acceptor_A"
                    if acceptor_g < len(true_path) and true_path[acceptor_g] == "Intron":
                        true_path[acceptor_g] = "Acceptor_G"

        for i in range(len(dna_seq)):
            nt = dna_seq[i]
            state = true_path[i]

            if state == "Skip":
                continue

            if state in states and nt in nt_order:
                nt_counts[state][nt] += 1

        for i in range(1, len(true_path)):
            prev_state = true_path[i - 1]
            cur_state = true_path[i]

            if prev_state == "Skip" or cur_state == "Skip":
                continue

            if prev_state in next_states:
                if cur_state in next_states[prev_state]:
                    state_change_counts[prev_state][cur_state] += 1

emission_log = {}
for state in states:
    emission_log[state] = {}
    total = 0
    for nt in nt_order:
        total += nt_counts[state][nt]

    for nt in nt_order:
        emission_log[state][nt] = math.log(nt_counts[state][nt] / total)

transition_log = {}
for state in states:
    transition_log[state] = {}
    total = 0
    for state2 in states:
        total += state_change_counts[state][state2]

    for state2 in states:
        if state_change_counts[state][state2] == 0:
            transition_log[state][state2] = -100
        else:
            transition_log[state][state2] = math.log(state_change_counts[state][state2] / total)

model_params = {
    "model_name": "splice_aware_HMM",
    "states": states,
    "transition_log": transition_log,
    "emission_log": emission_log
}

with open('model_params_splice.json', 'w') as f:
    json.dump(model_params, f, indent=4)
