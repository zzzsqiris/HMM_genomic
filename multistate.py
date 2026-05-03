import numpy
import math

states = ["Start", "Exon", "Intron", "End"]
n_s = len(states)

TP = [
    [0.1, 0.7, 0.1, 0.1],  # Start->Start/Exon/Intron/End
    [0.1, 0.7, 0.1, 0.1],  # Exon->Start/Exon/Intron/End
    [0.1, 0.1, 0.7, 0.1],  # Intron->Start/Exon/Intron/End
    [0.1, 0.1, 0.1, 0.7]   # End->Start/Exon/Intron/End
]

EP = [
    [0.25, 0.25, 0.25, 0.25], # Start: A/C/G/T
    [0.10, 0.40, 0.40, 0.10], # Exon: A/C/G/T
    [0.40, 0.10, 0.10, 0.40], # Intron: A/C/G/T
    [0.25, 0.25, 0.25, 0.25]  # End: A/C/G/T
]

# init by start
init = [0.7, 0.1, 0.1, 0.1]

# 0:A, 1:C, 2:G, 3:T
seq = [0, 2, 1, 2, 2, 3, 0, 0, 1, 2]

dpm = numpy.zeros((n_s, len(seq) + 1))
trace = numpy.full((n_s, len(seq) + 1), -1)

for s in range(n_s):
    dpm[s][0] = math.log(init[s])

# loop through seq
for i in range(1, len(seq)+1): 
    # loop through states
    for cur in range(n_s):
        # init as first state
        max_log_p = dpm[0][i-1] + math.log(TP[0][cur]) + math.log(EP[cur][seq[i-1]])
        best_prev_node = 0
        
        # loop through other states
        for prev in range(1, n_s):
            current_p = dpm[prev][i-1] + math.log(TP[prev][cur]) + math.log(EP[cur][seq[i-1]])
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


path = []
path.append(current)

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
