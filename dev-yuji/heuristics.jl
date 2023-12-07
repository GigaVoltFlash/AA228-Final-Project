include("../dev-emily/state.jl")

function heuristic_action(state, max_slew)
    # extract the state 
    koe = state.koe     
    attitude = state.attitude 
    dt  = state.dt            
    target_list = state.target_list  
    observed_list = state.observed_list    
    flag = .!Bool.(observed_list)  # flag = 1 if the target has not been observed
    id_targ = 1:length(observed_list)
    mu_E = 398600.4418  # [km^3/s^2]
    Re = 6375.0  # [km]

    for (idx, target_tup) in enumerate(target_list)
        if flag[idx] == 0
            continue
        end

        rv = koe2cart(copy(koe), mu_E)
        horizon_dist = sqrt( norm(rv[1:3])^2 - Re^2 )
        target_pos = ECEF_to_ECI([target_tup[1], target_tup[2], target_tup[3], 0, 0, 0], dt)
        # R_eci2rtn = ECI_to_RTN_matrix(rv)

        # look_vec_rtn = R_eci2rtn * (target_pos[1:3] .- rv[1:3]) 
        target_dis = norm(target_pos[1:3] .- rv[1:3]) 
        if horizon_dist < target_dis
            flag[idx] = 0
            continue
        end

        slews = get_slew_angle(koe, target_tup, dt) 
        if abs(attitude[1] - slews[1]) > max_slew || abs(attitude[2] - slews[2]) > max_slew
            flag[idx] = 0
        end

    end

    feas_target = target_list[flag.==1]
    feas_id_targ = id_targ[flag.==1]
    feas_target_ = [elem for tup in feas_target for elem in tup]
    feas_target_ = reshape(feas_target_, (length(feas_target), 6))
    # println("feas_target: ", size(feas_target))
    # println("feas_target: ", feas_target_[:,3])
    if length(feas_target_) == 0
        return 1
    end
    max_val, max_idx = findmax(feas_target_[:,3])
    # println("max_idx;", max_idx)
    max_id = feas_id_targ[max_idx]

    return max_id + 1 
end 

