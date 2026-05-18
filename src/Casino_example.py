import numpy
import math

states = ["F", "L"]
TP = [
    #[FF, FL]
    #[LF, LL]
    [0.9, 0.1],
    [0.2, 0.8]
]

EP = [
    [1/6, 1/6, 1/6,1 /6, 1/6, 1/6],     #F
    [0.1, 0.1, 0.1, 0.1, 0.1, 0.5]      #T
]

init = [0.5, 0.5]
seq = [0, 1, 3, 5, 5, 5, 5, 5, 5, 1, 4, 1, 1, 1,1, 1, 1, 1, 1, 1, 1]
dpm = numpy.zeros((2, len(seq)+1))
trace= numpy.full((2, len(seq)+1), -1)

dpm[0][0] = math.log(init[0])
dpm[1][0] = math.log(init[1])

for i in range (1, len(seq)+1):
    p00 = dpm[0][i-1] + math.log(TP[0][0]) + math.log(EP[0][seq[i-1]])
    p10 = dpm[1][i-1] + math.log(TP[1][0]) + math.log(EP[0][seq[i-1]])
    if p00 > p10:
        dpm[0][i] = p00
        trace[0][i] = 0
    else:
        dpm[0][i] = p10
        trace[0][i] = 1

    p01 = dpm[0][i-1] + math.log(TP[0][1]) + math.log(EP[1][seq[i-1]])
    p11 = dpm[1][i-1] + math.log(TP[1][1]) + math.log(EP[1][seq[i-1]])
    if p01 > p11:
        dpm[1][i] = p01
        trace[1][i] = 0
    else:
        dpm[1][i] = p11
        trace[1][i] = 1

print(dpm)
print(trace)

# get the path
# which row
if dpm[0][len(seq)] > dpm[1][len(seq)]:
    current_state = 0
else:
    current_state = 1
# loop back
path = []
path.append(current_state)
for i in range(len(seq), 0, -1):
    prev_state = trace[current_state][i]
    path.append(prev_state)
    current_state = prev_state

path.reverse()

path_names = []
for p in path:
    name = states[p]
    path_names.append(name)
print(path_names)