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

    for (idx, target_tup) in enumerate(target_list)
        if flag[idx] == 0
            continue
        end

        slews = get_slew_angle(koe, target_tup, dt)
        if abs(attitude[1] - slews[1]) < max_slew && abs(attitude[2] - slews[2]) < max_slew
            flag[idx] = 0
        end

    end

    feas_target = target_list[flag.==1]
    feas_id_targ = id_targ[flag.==1]
    feas_target_ = [elem for tup in feas_target for elem in tup]
    feas_target_ = reshape(feas_target_, (length(feas_target), 5))
    # println("feas_target: ", size(feas_target))
    # println("feas_target: ", feas_target_[:,3])
    max_val, max_idx = findmax(feas_target_[:,3])
    println("max_idx;", max_idx)
    max_id = feas_id_targ[max_idx]

    return max_id + 1 
end 

