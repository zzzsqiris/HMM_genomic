import HMM_utils
import os
import json
import math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("data_dir")
args = parser.parse_args()

data_dir = args.data_dir

# states
states = [
    "Exon",
    "Donor_1", "Donor_2", "Donor_3", "Donor_4", "Donor_5",
    "Intron",
    "Acceptor_1", "Acceptor_2", "Acceptor_3", "Acceptor_4", "Acceptor_5", "Acceptor_6"
]
nt_order = ['A', 'C', 'G', 'T']
fixed_nt = {
    "Donor_1": "G",
    "Donor_2": "T",
    "Donor_3": "A",
    "Donor_4": "A",
    "Donor_5": "G",
    "Acceptor_1": "T",
    "Acceptor_2": "T",
    "Acceptor_3": "T",
    "Acceptor_4": "C",
    "Acceptor_5": "A",
    "Acceptor_6": "G"
}

nt_counts = {}
state_change_counts = {}

# count tables
for state in states:
    nt_counts[state] = {'A': 1, 'C': 1, 'G': 1, 'T': 1}
    state_change_counts[state] = {}
    for state2 in states:
        state_change_counts[state][state2] = 0

# grammar
next_states = {
    "Exon": ["Exon", "Donor_1"],
    "Donor_1": ["Donor_2"],
    "Donor_2": ["Donor_3"],
    "Donor_3": ["Donor_4"],
    "Donor_4": ["Donor_5"],
    "Donor_5": ["Intron"],
    "Intron": ["Intron", "Acceptor_1"],
    "Acceptor_1": ["Acceptor_2"],
    "Acceptor_2": ["Acceptor_3"],
    "Acceptor_3": ["Acceptor_4"],
    "Acceptor_4": ["Acceptor_5"],
    "Acceptor_5": ["Acceptor_6"],
    "Acceptor_6": ["Exon"]
}

# transition start counts
for state in next_states:
    for state2 in next_states[state]:
        state_change_counts[state][state2] = 1

for filename in os.listdir(data_dir):
    if filename.endswith(".fa"):
        # input files
        fa_path = os.path.join(data_dir, filename)
        gff_path = os.path.join(data_dir, filename[:-3] + ".gff3")

        if not os.path.exists(gff_path):
            continue

        fasta_dict = HMM_utils.read_fasta(fa_path)
        gene_name = filename[:-3]
        dna_seq = fasta_dict[gene_name].upper()

        # true path
        true_path = HMM_utils.read_true_path(gff_path, len(dna_seq), False)
        features = HMM_utils.read_gff(gff_path)

        # splice states
        for ftype, start, end in features:
            if ftype == "intron":
                if end - start + 1 >= 11:
                    for j in range(5):
                        pos = start - 1 + j
                        if pos < len(true_path) and true_path[pos] == "Intron":
                            true_path[pos] = "Donor_" + str(j + 1)

                    for j in range(6):
                        pos = end - 6 + j
                        if pos < len(true_path) and true_path[pos] == "Intron":
                            true_path[pos] = "Acceptor_" + str(j + 1)

        # nt counts
        for i in range(len(dna_seq)):
            nt = dna_seq[i]
            state = true_path[i]

            if state == "Skip":
                continue

            if state in states and state not in fixed_nt and nt in nt_order:
                nt_counts[state][nt] += 1

        # transition counts
        for i in range(1, len(true_path)):
            prev_state = true_path[i - 1]
            cur_state = true_path[i]

            if prev_state == "Skip" or cur_state == "Skip":
                continue

            if prev_state in next_states:
                if cur_state in next_states[prev_state]:
                    state_change_counts[prev_state][cur_state] += 1

emission_prob = {}
for state in states:
    if state in fixed_nt:
        emission_prob[state] = {fixed_nt[state]: 1}
    else:
        emission_prob[state] = HMM_utils.count_to_prob(nt_counts[state])

transition_prob = {}
for state in states:
    transition_prob[state] = HMM_utils.count_to_prob(state_change_counts[state])

model_params = {
    "model_name": "long_splice_HMM",
    "states": states,
    "transition_prob": transition_prob,
    "emission_prob": emission_prob
}

with open('model_params_long_splice.json', 'w') as f:
    json.dump(model_params, f, indent=4)
