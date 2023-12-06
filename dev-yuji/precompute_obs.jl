using CSV, Random, DataFrames


df_traj = CSV.read("../src/traj_info.csv", DataFrame)
time = df_traj.t
x_ECI = df_traj.x_ECI
y_ECI = df_traj.y_ECI
z_ECI = df_traj.z_ECI

x_ECEF = df_traj.x_ECEF
y_ECEF = df_traj.y_ECEF
z_ECEF = df_traj.z_ECEF

mu_E = 398600.4418  # [km^3/s^2]
Re = 6375.0  # [km]

a_SC = 7057
horizon_dist = sqrt( a_SC^2 - Re^2 )

for (idx, target_tup) in enumerate(target_list)
    
    for (idx, t) in time
    
        target_pos = ECEF_to_ECI([target_tup[1], target_tup[2], target_tup[3], 0, 0, 0], t)
        target_dis = norm(target_pos[1:3] .- [x_ECI, y_ECI, z_ECI]) 

        slews = get_slew_angle(koe, target_tup, t) 
        if abs(attitude[1] - slews[1]) > max_slew || abs(attitude[2] - slews[2]) > max_slew
            flag[idx] = 0
        end
    
    end

end

