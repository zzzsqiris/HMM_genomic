import gzip
import math
import json

PAIRS = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "a": "t", "t": "a", "c": "g", "g": "c",
}

def complement (seq):
    comp_seq = ""
    if not isinstance(seq, str):
        seq = "".join(seq)
    for nt in seq:
        comp_seq += (PAIRS[nt])
    return comp_seq

def rev_comp(seq):
    if not isinstance(seq, str):
        seq = "".join(seq)
    comp_seq = complement(seq)
    rev_comp_seq = comp_seq[::-1]
    return rev_comp_seq

# open file
def smart_open(filepath):
    if filepath.endswith('.gz'):
        return gzip.open(filepath, mode='rt')
    else:
        return open(filepath, mode='rt')
    
# read fasta
def read_fasta(filepath):
    fasta_dict = {}
    seq_id = ""
    cur_seq = []

    with smart_open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    fasta_dict[seq_id] = "".join(cur_seq)
                seq_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if seq_id:
            fasta_dict[seq_id] = "".join(cur_seq)
            
    return fasta_dict

def read_gff(filepath):
    info = []
    with smart_open(filepath) as f:
        for line in f:
            parts = line.split('\t')
            source = parts[1]
            ftype = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            if source == "WormBase":
                info.append((ftype, start, end))
    return info

def read_pred_gff(pred_file):
    pred_info = []

    with smart_open(pred_file) as f:
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

# output label
def state_for_output(state):
    if state.startswith('Donor_') or state.startswith('Acceptor_'):
        return 'Intron'
    return state

# kmer
def get_kmer(seq, pos):
    if pos == 0:
        k = 1
    elif pos == 1:
        k = 2
    elif pos == 2:
        k = 3
    else:
        k = 4

    start = pos - k + 1
    return str(k), seq[start:pos + 1]

# true path
def read_true_path(gff_file, seq_length, is_multi_state):
    true_path = ['Intron'] * seq_length
    exon_label = [None] * seq_length
    intron_pos = [False] * seq_length
    gff_info = read_gff(gff_file)

    for ftype, start, end in gff_info:
        target_types = ['exon', 'five_prime_UTR', 'three_prime_UTR', 'CDS']
        if ftype in target_types:
            if is_multi_state:
                label = ftype
            else:
                label = 'Exon'
            for i in range(start - 1, end):
                if i < seq_length:
                    exon_label[i] = label
        elif ftype == 'intron':
            for i in range(start - 1, end):
                if i < seq_length:
                    intron_pos[i] = True

    for i in range(seq_length):
        if exon_label[i] is not None and intron_pos[i]:
            true_path[i] = 'Skip'
        elif exon_label[i] is not None:
            true_path[i] = exon_label[i]
        elif intron_pos[i]:
            true_path[i] = 'Intron'

    return true_path

# count to log
def count_to_log(counts):
    total = sum(counts.values())
    log_values = {}
    for nt, count in counts.items():
        if total > 0 and count > 0:
            prob = count / total
            log_values[nt] = math.log(prob)
        else:
            log_values[nt] = -100
    return log_values

# read model
def read_in_prob(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)
    all_states = data["states"]
    nt_order = ['A', 'C', 'G', 'T']

    if data.get("emission_type") == "kmer":
        EP = data["emission_log"]
    else:
        EP = []
        for state in all_states:
            state_ep = data["emission_log"][state]
            row = [state_ep[nt] for nt in nt_order]
            EP.append(row)

    TP = []
    for s_from in all_states:
        row = []
        for s_to in all_states:
            row.append(data["transition_log"][s_from][s_to])
        TP.append(row)
        
    return all_states, TP, EP
