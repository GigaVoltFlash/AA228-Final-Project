include("state.jl")

# Create an intial state
x0 = 0 # checl: are the test points centered around x=0?
y0 = 0
dydt = 1
alt = 50
attitude = (0, 0)
target_list, n_targets = create_target_list("src/obs_site.csv")

state = State(x0, y0, dydt, alt, attitude, target_list, zeros(n_targets))
# print(state)
# create action list
""" 
Action space:
1 -> Do nothing OR: rotate to neutral/some intermediate position
2 -> Image target 1
.
.
.
n+1 -> Image target n 
"""
A = 1:n_targets 

for t = 0:60
    global state
    state, reward = TR(state, rand(A))

end