import numpy
import math
import HMM_utils

states, TP, EP = HMM_utils.read_in_prob('model_params.json')
n_s = len(states)

fasta_dict = HMM_utils.read_fasta("smallgenes/ce.3.35.fa")
raw_seq = list(fasta_dict.values())[0].upper()
nt_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
seq = [nt_map[nt] for nt in raw_seq if nt in nt_map]


init_log = [math.log(0.5), math.log(0.5)]

dpm = numpy.zeros((n_s, len(seq) + 1))
trace = numpy.full((n_s, len(seq) + 1), -1)

for s in range(n_s):
    dpm[s][0] = init_log[s]

# loop through seq
for i in range(1, len(seq)+1): 
    obs = seq[i-1]
    # loop through states
    for cur in range(n_s):
        # init as first state
        max_log_p = dpm[0][i-1] + TP[0][cur] + EP[cur][obs]
        best_prev_node = 0
        
        # loop through other states
        for prev in range(1, n_s):
            current_p = dpm[prev][i-1] + TP[prev][cur] + EP[cur][obs]
            if current_p > max_log_p:
                max_log_p = current_p
                best_prev_node = prev
        
        dpm[cur][i] = max_log_p
        trace[cur][i] = best_prev_node

# find last state, init as first state
max_val = dpm[0][len(seq)]
current = 0

for s in range(1, n_s):
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
for p in path:
    name = states[p]
    path_names.append(name)
print(path_names)
