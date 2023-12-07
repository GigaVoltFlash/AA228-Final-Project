println("Getting required packages")
using Plots
using CSV, Random, DataFrames
using GLMakie
using GeoDatasets
using JLD2
include("../dev-yuji/heuristics.jl")
include("../src/state.jl")

MARKER_SIZE_FACTOR = 25 # CHANGE TO WHAT LOOKS GOOD
# constants
Re = 6371;
mu_E = 3.986004418e5; #km^3/m^2  
slew_limit = 15

println("Loading in data")
@load "heuristic_runs.jld2" state_list action_list cross_ang_des along_ang_des observed_target_lambda observed_target_phi observed_target_reward num_observed reward_totals target_lambda target_phi target_mean_reward
@load "mcts_random_runs.jld2" state_list_mcts action_list_mcts cross_ang_des_mcts along_ang_des_mcts mcts_observed_target_lambda mcts_observed_target_phi mcts_observed_target_reward mcts_num_observed mcts_reward_totals target_lambda target_phi target_mean_reward

println("Doing some last minute angle calcs")
# Get some angles stuff
x_angle_list = []
y_angle_list = []
index = 1:length(state_list)
for i in 1:length(state_list)
    push!(x_angle_list, state_list[i].attitude[1])
    push!(y_angle_list, state_list[i].attitude[2])
end

x_angle_list_mcts = []
y_angle_list_mcts = []
index_mcts = 1:length(state_list_mcts)
for i in 1:length(state_list_mcts)
    push!(x_angle_list_mcts, state_list_mcts[i].attitude[1])
    push!(y_angle_list_mcts, state_list_mcts[i].attitude[2])
end

println("Plotting")
# Plots.plot(index, action_list, label="Heuristic")
# Plots.plot!(index, action_list_mcts, title="Actions with MCTS", label="MCTS", show=true)
# savefig("action_plot.png")

# random_along = Plots.plot(index, [y_angle_list, along_ang_des], label=["Actual Angle" "Desired Angle"], title="Along track angle for Baseline")
# mcts_along = Plots.plot(index, [y_angle_list_mcts, along_ang_des_mcts], label=["Actual Angle" "Desired Angle"], title="Along track angle for MCTS")
# Plots.plot(random_along, mcts_along,layout=(1, 2), size=(1000, 400), ylabel="Degrees", xlabel="Time (s)", show=true)
# savefig("along_angles_plot.png")

# random_cross = Plots.plot(index, [x_angle_list, cross_ang_des], label=["State Heuristic" "Desired Heuristic"])
# mcts_cross = Plots.plot(index, [x_angle_list_mcts, cross_ang_des_mcts], label=["State MCTS" "Desired MCTS"])
# Plots.plot(random_cross, mcts_cross,layout=(1, 2), title="Cross track angle", size=(800, 400), show=true)
# savefig("cross_angles_plot.png")


# ######### GLOBE PLOTTER #########
koe = [state.koe for state in state_list]
dts = [state.dt for state in state_list]
state_ECEF = zeros(6, length(dts))
state_ECI = zeros(6, length(dts))
state_lambda = zeros(length(dts))
state_phi = zeros(length(dts))
for (i, koe) in enumerate(koe)
    rv = koe2cart(koe, mu_E)
    rv_ecef = ECI_to_ECEF(rv, dts[i])
    state_lambda[i] = atand(rv_ecef[2], rv_ecef[1])
    state_phi[i] = atand(rv_ecef[3], sqrt(rv_ecef[1]^2 + rv_ecef[2]^2))
end

layout = @layout [a b]
f = Figure(size = (1000, 800), layout=layout)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
ax1 = Axis(f[1, 1], aspect=DataAspect(), title = "Heuristic")
sca1 =GLMakie.scatter!(ax1, target_lambda, target_phi, overdraw = true, color = :grey, markersize = MARKER_SIZE_FACTOR*target_mean_reward)
sca2 =GLMakie.scatter!(ax1, observed_target_lambda , observed_target_phi, overdraw = true, color = :blue, markersize = MARKER_SIZE_FACTOR*observed_target_reward)
GLMakie.contour!(ax1, lon, lat, data, levels=[0.5], color=:black, linewidth=0.5)
path = lines!(ax1, state_lambda, state_phi, overdraw = true, color = :blue, linewidth = 2,)
GLMakie.ylims!(ax1, low = -66, high = 66)
GLMakie.xlims!(ax1, low = -140, high = -85)

mcts_koe = [state.koe for state in state_list_mcts]
mcts_dts = [state.dt for state in state_list_mcts]
state_ECEF = zeros(6, length(mcts_dts))
state_ECI = zeros(6, length(mcts_dts))
mcts_state_lambda = zeros(length(mcts_dts))
mcts_state_phi = zeros(length(mcts_dts))
for (i, koe) in enumerate(mcts_koe)
    rv = koe2cart(koe, mu_E)
    rv_ecef = ECI_to_ECEF(rv, mcts_dts[i])
    mcts_state_lambda[i] = atand(rv_ecef[2], rv_ecef[1])
    mcts_state_phi[i] = atand(rv_ecef[3], sqrt(rv_ecef[1]^2 + rv_ecef[2]^2))
end


k = 10.580 
Δλ = 10
p_br = [-110 + k + Δλ, -60]
p_bl = [-110 + k - Δλ, -60]
p_tr = [-110 - k + Δλ,  60]
p_tl = [-110 - k - Δλ,  60]

l1 = lines!(ax1, [p_br[1], p_bl[1]], [p_br[2], p_bl[2]], overdraw=true, color=:orange, linestyle=:solid, linewidth=2,)
l2 = lines!(ax1, [p_tr[1], p_tl[1]], [p_tr[2], p_tl[2]], overdraw=true, color=:orange, linestyle=:solid, linewidth=2,)
l3 = lines!(ax1, [p_br[1], p_tr[1]], [p_br[2], p_tr[2]], overdraw=true, color=:orange, linestyle=:solid, linewidth=2,)
l4 = lines!(ax1, [p_bl[1], p_tl[1]], [p_bl[2], p_tl[2]], overdraw=true, color=:orange, linestyle=:solid, linewidth=2,)


lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
ax2 = Axis(f[1, 2], aspect=DataAspect(), title = "MCTS with Random Rollouts")
GLMakie.scatter!(ax2, target_lambda, target_phi, overdraw = true, color = :grey, markersize = MARKER_SIZE_FACTOR*target_mean_reward)
GLMakie.scatter!(ax2, mcts_observed_target_lambda , mcts_observed_target_phi, overdraw = true, color = :blue, markersize = MARKER_SIZE_FACTOR*mcts_observed_target_reward)
GLMakie.contour!(ax2, lon, lat, data, levels=[0.5], color=:black, linewidth=0.5)
lines!(ax2, mcts_state_lambda, mcts_state_phi, overdraw = true, color = :blue, linewidth = 2,)
GLMakie.ylims!(ax2, low = -66, high = 66)
GLMakie.xlims!(ax2, low = -140, high = -85)

l1 = lines!(ax2, [p_br[1], p_bl[1]], [p_br[2], p_bl[2]], overdraw=true, color=:orange, linestyle=:solid, linewidth=2,)
l2 = lines!(ax2, [p_tr[1], p_tl[1]], [p_tr[2], p_tl[2]], overdraw=true, color=:orange, linestyle=:solid, linewidth=2,)
l3 = lines!(ax2, [p_br[1], p_tr[1]], [p_br[2], p_tr[2]], overdraw=true, color=:orange, linestyle=:solid, linewidth=2,)
l4 = lines!(ax2, [p_bl[1], p_tl[1]], [p_bl[2], p_tl[2]], overdraw=true, color=:orange, linestyle=:solid, linewidth=2,)

legend = Legend(f[1,3], [sca1, sca2, path, l1], ["Unobserved targets", "Observed targets", "Satellite Path", "Target Bounding Box"], framevisible=true)
GLMakie.save("./map_plot.png", f)

