import random

def get_next_state(current_state):
    i = random.random()
    if current_state == "F":
        if i <= 0.8:
            return "F"
        else:
            return "L"
    elif current_state == "L":
        if i <= 0.9:
            return "L"
        else:
            return "F"

def emission(state):
    if state == "F":
        return random.randint(1, 6)
    elif state == "L":
        i = random.randint(1, 10)
        if i <= 5:
            return i
        else:
            return 6

current_state = random.choice(["F", "L"])

state_record = []
obs_record = []
for _ in range(100):
    current_state = get_next_state(current_state)
    obs = emission(current_state)
    state_record.append(current_state)
    obs_record.append(obs)

print(state_record)
print(obs_record)