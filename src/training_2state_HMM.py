import HMM_utils
import os
import json
import math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("data_dir")
args = parser.parse_args()

data_dir = args.data_dir
exon_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
intron_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
total_exon_len = 0
num_exons = 0
total_intron_len = 0
num_introns = 0

def add_segment(state, length):
    global total_exon_len, num_exons, total_intron_len, num_introns #### no global

    if state == "Exon":
        total_exon_len += length
        num_exons += 1
    elif state == "Intron":
        total_intron_len += length
        num_introns += 1

for filename in os.listdir(data_dir):
    # read gff3 and fasta file
    if filename.endswith(".fa"):            ### if not endswith, less indent
        fa_path = os.path.join(data_dir, filename)
        fasta_dict = HMM_utils.read_fasta(fa_path)
        gff_path = os.path.join(data_dir, filename[:-3] + ".gff3")
        
        gene_name = filename[:-3]
        dna_seq = fasta_dict[gene_name]
        true_path = HMM_utils.read_true_path(gff_path, len(dna_seq), False)

        # counting for EP
        for i in range(len(dna_seq)):
            nt = dna_seq[i].upper()
            state = true_path[i]

            if state == 'Exon':
                if nt in exon_counts:
                    exon_counts[nt] += 1
            elif state == 'Intron':
                if nt in intron_counts:
                    intron_counts[nt] += 1

        # counting average state length for TP
        current_state = ""
        current_len = 0

        for state in true_path:
            if state == "Skip":
                if current_len > 0:
                    add_segment(current_state, current_len)
                current_state = ""
                current_len = 0
            elif state == current_state:
                current_len += 1
            else:
                if current_len > 0:
                    add_segment(current_state, current_len)
                current_state = state
                current_len = 1

        if current_len > 0:
            add_segment(current_state, current_len)

ep_exon = HMM_utils.count_to_prob(exon_counts)
ep_intron = HMM_utils.count_to_prob(intron_counts)

avg_exon_len = total_exon_len / num_exons
p_ei = 1 / avg_exon_len
p_ee = 1 - p_ei

avg_intron_len = total_intron_len / num_introns
p_ie = 1 / avg_intron_len
p_ii = 1 - p_ie

model_params = {
    "model_name": "E_I_two_state_HMM",
    "states": ["Exon", "Intron"],
    "transition_prob": {
        "Exon": {"Exon": p_ee, "Intron": p_ei},
        "Intron": {"Exon": p_ie, "Intron": p_ii}
    },
    "emission_prob": {
        "Exon": ep_exon,
        "Intron": ep_intron
    }
}

with open('model_params_2state.json', 'w') as f:
    json.dump(model_params, f, indent=4)
