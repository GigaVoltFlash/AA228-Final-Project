using CSV, Random, DataFrames
include("../src/state.jl")


df_traj = CSV.read("src/traj_info.csv", DataFrame)
df_target = CSV.read("src/obs_site_Earth_50.csv", DataFrame)

time = df_traj.t
n_targets = length(df_target.id)
x_ECI = df_traj.x_ECI
y_ECI = df_traj.y_ECI
z_ECI = df_traj.z_ECI
vx_ECI = df_traj.vx_ECI
vy_ECI = df_traj.vy_ECI
vz_ECI = df_traj.vz_ECI

x_ECEF = df_traj.x_ECEF
y_ECEF = df_traj.y_ECEF
z_ECEF = df_traj.z_ECEF

target_pos = cat(df_target.x, df_target.y, df_target.z, dims = 2)

mu_E = 398600.4418  # [km^3/s^2]
Re = 6375.0  # [km]

a_SC = 7057
horizon_dist = sqrt( a_SC^2 - Re^2 )
println("horizon_dist: ", horizon_dist)

dt = 30 
t_min = zeros(n_targets)  
t_max = ones(n_targets) * 2000
max_slew = 10

for (idx, t) in enumerate(time)
    
    target_pos = Matrix(df_target[:, [:x, :y, :z]])
    r_ECEF = [x_ECEF[idx], y_ECEF[idx], z_ECEF[idx]]
    # println("r_ECEF: ", r_ECEF)
    # println("target_pos: ", target_pos)
    a = target_pos .- transpose(r_ECEF)
    # println("a: ", a)
    target_dist = sqrt.(a[:,1].^2 + a[:,1].^2 + a[:,1].^2)
    # println("target_dist: ", target_dist)
    target_dis = sqrt.(a[:,1].^2 + a[:,1].^2 + a[:,1].^2) .< horizon_dist
    # println("target_dis: ", target_dis)

    for target_tup in eachrow(df_target)
        # println("target_tup: ", target_tup)
        id = target_tup.id
        # slew = w.r.t. the nadir direction (i.e., RTN)
        slews = get_slew_angle2([x_ECI[idx], y_ECI[idx], z_ECI[idx], vx_ECI[idx], vy_ECI[idx], vz_ECI[idx]], [target_tup.x, target_tup.y, target_tup.z], t) 
        slew_cond = abs(slews[1]) < max_slew && abs(slews[2]) < max_slew
        dist_cond = target_dis[id]
        observable = slew_cond && dist_cond
        
        if observable && t_min[id] == 0 
            t_min[id] = t
        elseif !observable && t_min[id] > 0
            t_max[id] = t - dt 
        end  
    
    end

end

t_max[t_min .== 0 .&& t_max.==2000] .= 0

df_target.t_min = t_min
df_target.t_max = t_max

CSV.write("target_info.csv", df_target)

