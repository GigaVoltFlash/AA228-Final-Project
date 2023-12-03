include("../dev-emily/state.jl")

function heuristic_action(state, max_slew)
    # extract the state 
    koe = state.koe     
    attitude = state.attitude 
    dt  = state.dt            
    target_list = state.target_list  
    observed_list = state.observed_list    
    flag = .!observed_list  # flag = 1 if the target has not been observed
    id_targ = 1:length(observed_list)

    for (idx, target_tup) in enumerate(target_list)
        if flag[idx] == 0
            continue
        end

        slews = get_slew_angle(koe, target_tup, dt_JD)
        if abs(atttiude[1] - slews[1]) < max_slew && abs(attitude[2] - slews[2]) < max_slew
            flag[i] = 0
        end

    end

    feas_target  = target_list[flag.==1]
    feas_id_targ = id_targ[flag.==1]
    max_val, max_idx = findmax(feas_target[:, 3])
    max_id = feas_id_targ[max_idx]

    return max_id + 1 
end 

