println("Getting required packages")
using Plots
using CSV, Random, DataFrames
using GLMakie
using GeoDatasets
using JLD2

include("../src/state.jl")

println("Setting up problem")
# SET THE DATASET YOU WANT TO TEST
dataset = "../src/obs_site_Earth_1000.csv"
# SET THE NUMBER OF TIMES YOU WANT TO RUN THIS
num_runs = 1

target_list, n_targets = create_target_list_3d(dataset)


# constants
Re = 6371;
mu_E = 3.986004418e5; #km^3/m^2  
slew_limit = 15

# MCTS parameters
gamma = 0.95 # discount factor for the MDP problem
m = 500 # simulation count
c = 4000 # exploration constant
d = 30 # simulation depth
d_r = 15 # rollout depth
time_step = 30

# reading in dataset for plotting in the future
df = CSV.read(dataset, DataFrame)
target_lambda = df.lambda
target_phi = df.phi
target_mean_reward = df.r_mean


####################################### MCTS functions #######################################
# A shortened representation of MDP that works for online planning (no full state-space, TR captures T and R)
struct MDP 
    gamma # Discount factor
    A # Action space
    TR # Transition model
end 

function rollout(P, s, d)
    s_ = copy(s)
    ret = 0.0 
    for t in 1:d
        a = rand(P.A) # FOR RANDOM MCTS
        s_, r = P.TR(s_, a)
        ret += P.gamma^(t-1) * r
    end 
    return ret 
end

struct MonteCarloTreeSearch
    P
    N   # visit count
    Q   # action-value estimate
    d   # depth to go through for the Monte Carlo search
    d_r # depth to go through the rollout at the end of the MCTS
    m   # num of simulations
    c   # exploration constant 
    U   # value function estimate 
end 

function (π::MonteCarloTreeSearch)(s)
    for k in 1:π.m
        simulate!(π,s)
    end 
    return argmax(a->π.Q[(s,a)], π.P.A)
end 

function simulate!(π::MonteCarloTreeSearch, s, d=π.d)
    if d <= 0 
        return rollout(π.P, s, π.d_r)
    end 
    P, N, Q, c = π.P, π.N, π.Q, π.c
    A, TR, gamma = P.A, P.TR, P.gamma
    if !haskey(N,(s,first(A)))  # if (s,a) has never been visited
        for a in A
            N[(s,a)] = 0 
            Q[(s,a)] = 0
        end 
        return rollout(π.P, s, π.d_r)   
    end
    
    # if (s,a) has been visited
    a = explore(π,s)
    # println("State before TR in simulate \n", s)
    s_, r = TR(s,a)
    # println("State after TR in simulate \n", s)
    q = r + gamma* simulate!(π, s_, d-1) 
    N[(s,a)] += 1  # +1 for visit count 
    Q[(s,a)] += (q-Q[(s,a)]) / N[(s,a)]  # The more you visited, the update of the Q will (usually) converge 
    return q 
end 

bonus(Nsa, Ns) = Nsa==0 ? Inf : sqrt(log(Ns)/Nsa)

function explore(π::MonteCarloTreeSearch, s)
    A, N, Q, c = π.P.A, π.N, π.Q, π.c
    Ns = sum(N[(s,a)] for a in A)
    # objective = Q+bonus term 
    # if there is no past visit, then that exploration is always prioritized 
    return argmax(a -> Q[(s,a)] + c*bonus(N[(s,a)], Ns), A)  
end 

####################### RUNNING MCTS ####################### 
mcts_reward_totals = zeros(num_runs)
mcts_num_observed = zeros(num_runs)
action_space =  1:n_targets+1 

P = MDP(gamma, action_space, TR_orbit)

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

state_list_mcts = []
action_list_mcts = []
# Create an intial state
observed_list = zeros(n_targets)
state = State3d(koe0, att, dt, target_list, observed_list)

for j in 1:num_runs
    global state_list_mcts, action_list_mcts, state
    observed_list = zeros(n_targets)
    state = State3d(koe0, att, dt, target_list, observed_list)

    Q = Dict{Tuple{State3d, Int}, Float64}()
    N = Dict{Tuple{State3d, Int}, Int}()
    U = Dict{State, Float64}()

    # Q, N, and U are going to be difficult to define
    mcts_run = MonteCarloTreeSearch(P, N, Q, d, d_r, m, c, U)

    state_list_mcts = []
    action_list_mcts = []
    for t = 0:time_step:2000
        global state_list_mcts, action_list_mcts, state
        # Push to the full list
        state_copy = copy(state)
        push!(state_list_mcts, state_copy)
        # println("State copy")
        # println(state_copy.observed_list)

        # Run MCTS
        state_for_mcts = copy(state)
        mcts_action = mcts_run(state_for_mcts)

        # Use the action
        push!(action_list_mcts, mcts_action)
        state, reward = TR_orbit(state, mcts_action, time_step)
        mcts_reward_totals[j] += reward

        # Calculate progress percentage
        progress = t / 2000 * 100

        # Create a string representation of the progress bar
        bar_length = 20
        bar_progress = round(Int, progress / (100 / bar_length))
        bar_string = "[" * repeat("=" , bar_progress) * repeat(" " , bar_length - bar_progress) * "]"

        # Print the progress bar
        print("\rMCTS with Random Rollouts Progress $j: $bar_string $(round(progress, digits=2))%")

        # Sleep for a short duration to see the progress bar
        sleep(0.1)
    end
    mcts_num_observed[j] = sum(state.observed_list)
    println("")
end

println("##################### RESULTS #####################")
println("Mean reward for MCTS with random rollouts : ", mean(mcts_reward_totals))
println("Mean no. of observations for MCTS with random rollouts : ", mean(mcts_num_observed))
println("##################### RESULTS #####################")

println("Gathering data for plotting")
# Gathering data for plotting
mcts_observed_target_lambda = []
mcts_observed_target_phi = []
mcts_observed_target_reward = []
for i in 1:length(state.observed_list)
    if state.observed_list[i] == 1
        push!(mcts_observed_target_lambda, target_lambda[i])
        push!(mcts_observed_target_phi, target_phi[i])
        push!(mcts_observed_target_reward, target_mean_reward[i])
    end
end
mcts_observed_target_lambda = Float64.(mcts_observed_target_lambda)
mcts_observed_target_phi = Float64.(mcts_observed_target_phi)
mcts_observed_target_reward = Float64.(mcts_observed_target_reward)

angs_list_mcts = []
for i in 1:length(action_list_mcts)
    # println(i)
    action = action_list_mcts[i]
    # println(action)
    state = state_list_mcts[i]
    # println(state)
    if action > 1
        target = state.target_list[action-1]
        push!(angs_list_mcts, get_slew_angle(state.koe, target, state.dt))
    else 
        if i > 1
            push!(angs_list_mcts, angs_list_mcts[i-1])
        else
            push!(angs_list_mcts, (0,0))
        end
    end
end
angs_list_mcts = collect(angs_list_mcts)
cross_ang_des_mcts = [t[1] for t in angs_list_mcts]
along_ang_des_mcts = [t[2] for t in angs_list_mcts];

println("Saving variables to mcts_random_runs.jld2")
@save "mcts_random_runs1000.jld2" state_list_mcts action_list_mcts cross_ang_des_mcts along_ang_des_mcts mcts_observed_target_lambda mcts_observed_target_phi mcts_observed_target_reward mcts_num_observed mcts_reward_totals target_lambda target_phi target_mean_reward
println("DONE")