println("Getting required packages")
using Plots
using CSV, Random, DataFrames
using GLMakie
using GeoDatasets
using JLD2

include("../dev-yuji/heuristics.jl")
include("state.jl")

println("Setting up problem")
# SET THE DATASET YOU WANT TO TEST
dataset = "obs_site_Earth_200.csv"
# SET THE NUMBER OF TIMES YOU WANT TO RUN THIS
num_runs = 10

# constants
Re = 6371;
mu_E = 3.986004418e5; #km^3/m^2  
slew_limit = 15


# reading in dataset for plotting in the future
df = CSV.read(dataset, DataFrame)
target_lambda = df.lambda
target_phi = df.phi
target_mean_reward = df.r_mean

################ RUNNING GREEDY HEURISTIC ####################
println("Running greedy heuristic")

# setup parameters
reward_totals = zeros(num_runs)
num_observed = zeros(num_runs)
target_list, n_targets = create_target_list_3d(dataset)
time_step = 30 # seconds
num_minutes = 30
slew_limit = 15



# Orbit propogation params
koe0 = [
    7057.0,  # a [km]
    0.000879,  # e
    deg2rad(98.12),  # i
    deg2rad(255.09),  # Ω [rad]
    deg2rad(225.37),  # ω [rad]
    deg2rad(75.0),  # M [rad]
]
att = [0.,0.]
dt = 0

observed_list = zeros(n_targets)
state = State3d(koe0, att, dt, target_list, observed_list)
state_list = []
action_list = []
for j in 1:num_runs
    global state, state_list, action_list
    # Create initial state
    A = 1:n_targets+1 
    observed_list = zeros(n_targets)
    state = State3d(koe0, att, dt, target_list, observed_list)
    state_list = []
    action_list = []
    for t in 0:time_step:2000
        global state, state_list, action_list
        state_copy = copy(state)
        push!(state_list, state_copy)
        action = heuristic_action(state, slew_limit)
        push!(action_list, action)
        state, reward = TR_orbit(state, action, time_step)
        reward_totals[j] += reward

        # Calculate progress percentage
        progress = t / 2000 * 100

        # Create a string representation of the progress bar
        bar_length = 20
        bar_progress = round(Int, progress / (100 / bar_length))
        bar_string = "[" * repeat("=" , bar_progress) * repeat(" " , bar_length - bar_progress) * "]"

        # Print the progress bar
        print("\rBaseline Heuristic Progress $j: $bar_string $(round(progress, digits=2))%")

        # Sleep for a short duration to see the progress bar
        sleep(0.1)
    end
    num_observed[j] = sum(state.observed_list)
    println("")
end

println("##################### RESULTS #####################")
println("Final reward for baseline heuristic: ", mean(reward_totals))
println("Median number of sites observed for baseline heuristic ", median(num_observed))
println("##################### RESULTS #####################")

println("Gathering data for plotting")

# Gathering data for heuristic from the last state
observed_target_lambda = []
observed_target_phi = []
observed_target_reward = []
for i in 1:length(state.observed_list)
    if state.observed_list[i] == 1
        push!(observed_target_lambda, target_lambda[i])
        push!(observed_target_phi, target_phi[i])
        push!(observed_target_reward, target_mean_reward[i])
    end
end
observed_target_lambda = Float64.(observed_target_lambda)
observed_target_phi = Float64.(observed_target_phi)
observed_target_reward = Float64.(observed_target_reward)

angs_list = []
for i in 1:length(action_list)
    # println(i)
    action = action_list[i]
    # println(action)
    state = state_list[i]
    if action > 1
        target = state.target_list[action-1]
        push!(angs_list, get_slew_angle(state.koe, target, state.dt))
    else 
        if i > 1
            push!(angs_list, angs_list[i-1])
        else
            push!(angs_list, (0,0))
        end
    end
end
angs_list = collect(angs_list)
cross_ang_des = [t[1] for t in angs_list];
along_ang_des = [t[2] for t in angs_list];

println("Saving variables to heuristic_runs.jld2")
@save "heuristic_runs.jld2" state_list action_list cross_ang_des along_ang_des observed_target_lambda observed_target_phi observed_target_reward num_observed reward_totals target_lambda target_phi target_mean_reward

println("DONE")