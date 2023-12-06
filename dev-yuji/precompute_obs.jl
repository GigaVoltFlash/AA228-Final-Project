using CSV, Random, DataFrames
include("../dev-emily/state.jl")


df_traj = CSV.read("../src/traj_info.csv", DataFrame)
df_target = CSV.read("../src/obs_site_Earth2.csv", DataFrame)

time = df_traj.t
x_ECI = df_traj.x_ECI
y_ECI = df_traj.y_ECI
z_ECI = df_traj.z_ECI

x_ECEF = df_traj.x_ECEF
y_ECEF = df_traj.y_ECEF
z_ECEF = df_traj.z_ECEF

target_pos = df_target[:,[:x, :y, :z]]

mu_E = 398600.4418  # [km^3/s^2]
Re = 6375.0  # [km]

a_SC = 7057
horizon_dist = sqrt( a_SC^2 - Re^2 )

dt = 30 
t_min = 0 
t_max = 2000
max_slew = 15

for (idx, t) in time
    
    r_ECEF = [x_ECEF[idx], y_ECEF[idx], z_ECEF[idx]]
    target_dis = norm(target_pos .- r_ECEF) .< horizon_dist

    for (j, target_tup) in enumerate(target_list)
        
        # id = target_tup.id
        
        slews = get_slew_angle2([x_ECI[idx], y_ECI[idx], z_ECI[idx]], target_tup, t) 
        slew_cond = abs(attitude[1] - slews[1]) < max_slew && abs(attitude[2] - slews[2]) < max_slew
        dist_cond = target_dis[j]


        
        if (abs(attitude[1] - slews[1]) > max_slew || abs(attitude[2] - slews[2]) > max_slew) && t_min > 0 
            t_max = t_min - dt 
        elseif t_min == 0 && horizon_dist < target_dis
            t_min = t 
        end        


    
    end

end

