import HMM_utils
import os
import json
import math

data_dir = "smallgenes"
exon_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
intron_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
total_exon_len = 0
num_exons = 0
total_intron_len = 0
num_introns = 0

for filename in os.listdir(data_dir):
    # read gff3 and fasta file
    if filename.endswith(".fa"):
        fa_path = os.path.join(data_dir, filename)
        fasta_dict = HMM_utils.read_fasta(fa_path)
        gff_path = os.path.join(data_dir, filename[:-3] + ".gff3")
        
        gene_name = filename[:-3]
        dna_seq = fasta_dict[gene_name]
        features = HMM_utils.read_gff(gff_path)

        # counting for EP
        for ftype, start, end in features:
            sub_seq = dna_seq[start-1 : end].upper()
            if ftype == 'exon':
                for nt in sub_seq:
                    if nt in exon_counts:
                        exon_counts[nt] += 1
                total_exon_len += (end - start + 1)
                num_exons += 1
            elif ftype == 'intron':
                for nt in sub_seq:
                    if nt in intron_counts:
                        intron_counts[nt] += 1
                total_intron_len += (end - start + 1)
                num_introns += 1

ep_exon = HMM_utils.log_probs(exon_counts)
ep_intron = HMM_utils.log_probs(intron_counts)

avg_exon_len = total_exon_len / num_exons
p_ei = 1 / avg_exon_len
p_ee = 1 - p_ei

avg_intron_len = total_intron_len / num_introns
p_ie = 1 / avg_intron_len
p_ii = 1 - p_ie

model_params = {
    "model_name": "E_I_two_state_HMM",
    "states": ["Exon", "Intron"],
    "transition_log": {
        "Exon": {"Exon": math.log(p_ee), "Intron": math.log(p_ei)},
        "Intron": {"Exon": math.log(p_ie), "Intron": math.log(p_ii)}
    },
    "emission_log": {
        "Exon": ep_exon,
        "Intron": ep_intron
    }
}

with open('model_params.json', 'w') as f:
    json.dump(model_params, f, indent=4)