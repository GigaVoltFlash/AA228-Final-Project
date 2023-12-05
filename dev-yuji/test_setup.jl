include("../dev-emily/state.jl")
include("heuristics.jl")

# create initial 3d state
# koe = [6571, 0.001, 0, 0, 0, 0]

koe = [
    7057.0,  # a [km]
    0.000879,  # e 
    deg2rad(98.12),  # i
    deg2rad(255.09),  # Ω [rad]
    deg2rad(225.37),  # ω [rad]
    deg2rad(75.0),  # M [rad]
]
att = [0.0,0.0]
dt = 0
target_list, n_targets = create_target_list_3d("src/obs_site_Earth.csv")
observed_list = zeros(n_targets)
state = State3d(koe, att, dt, target_list, observed_list)

""" 
Action space:
1 -> Do nothing OR: rotate to neutral/some intermediate position
2 -> Image target 1
.
.
.
n+1 -> Image target n 
"""
A = 1:n_targets + 1

time_step = 30 # seconds
reward_total = 0
for t in 0:time_step:30*60
    global state, reward_total
    # a = rand(A)
    a = heuristic_action(state, 15)
    print("Choosing action ")
    print(a)
    print("\n")

    state, reward = TR_orbit(state, a, time_step)
    reward_total += reward

end

print("Total reward: ")
println(reward_total)
