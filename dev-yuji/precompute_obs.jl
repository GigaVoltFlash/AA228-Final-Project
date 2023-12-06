using CSV, Random, DataFrames
include("../dev-emily/state.jl")


df_traj = CSV.read("src/traj_info.csv", DataFrame)
df_target = CSV.read("src/obs_site_Earth_50.csv", DataFrame)

time = df_traj.t
n_targets = length(df_target.id)
println("t:", time)
x_ECI = df_traj.x_ECI
y_ECI = df_traj.y_ECI
z_ECI = df_traj.z_ECI

x_ECEF = df_traj.x_ECEF
y_ECEF = df_traj.y_ECEF
z_ECEF = df_traj.z_ECEF

target_pos = cat(df_target.x, df_target.y, df_target.z, dims = 2)

mu_E = 398600.4418  # [km^3/s^2]
Re = 6375.0  # [km]

a_SC = 7057
horizon_dist = sqrt( a_SC^2 - Re^2 )

dt = 30 
t_min = zeros(length(time))  
t_max = ones(length(time)) * 2000
max_slew = 15

for (idx, t) in enumerate(time)
    
    target_pos = Matrix(df_target[:, [:x, :y, :z]])
    r_ECEF = [x_ECEF[idx], y_ECEF[idx], z_ECEF[idx]]
    a = target_pos .- transpose(r_ECEF)
    target_dis = sqrt.(a[:,1].^2 + a[:,1].^2 + a[:,1].^2) .< horizon_dist

    for (j, target_tup) in enumerate(target_list)
        
        id = target_tup.id
        slews = get_slew_angle2([x_ECI[idx], y_ECI[idx], z_ECI[idx]], target_tup, t) 
        slew_cond = abs(attitude[1] - slews[1]) < max_slew && abs(attitude[2] - slews[2]) < max_slew
        dist_cond = target_dis[j]
        observable = slew_cond && dist_cond
        
        if observable && t_min[id] == 0 
            t_min[id] = t
        elseif !observable && t_min[id] > 0
            t_max[id] = t - dt 
        end  
    
    end

end

df_target.t_min = t_min
df_target.t_max = t_max

CSV.write("target_info.csv", df_target)

