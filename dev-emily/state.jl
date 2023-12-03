using DataFrames
using CSV
using Distributions
using LinearAlgebra

include("../src/orbits.jl")

mutable struct State
    # lat             # in degrees
    # lon             # in degrees
    x               # for the simple test case only - in km
    y               # for the simple test case only - in km
    dydt            # "velocity" of the s/c -- how far does it move in a time step. Probably can be defined outside this struct
    alt             # in km(?)
    attitude        # Tuple of angles - for now: (x_angle, y_angle)
    target_list     # n entries in list; each entry is a tuple of (x, y, r)
    observed_list   # n entries in list; 1 indicates if a target was observed
end

Base.copy(state::State) = State(deepcopy(state.x), deepcopy(state.y), deepcopy(state.dydt), deepcopy(state.alt), deepcopy(state.attitude), deepcopy(state.target_list), deepcopy(state.observed_list))

mutable struct State3d
    koe             # [a,e,i,Ω,ω,M] w/ angles in radians, distance in km 
    attitude        # vector of angles - for now: [cross_track_angle (+ to left), along_track_angle (+ in vel direction)].
    dt              # time elapsed since beginning of scenario
    target_list     # n entries in list; each entry is a tuple of (x, y, z, r_mean, r_cov) w/ positions in ECEF
    observed_list   # n entries in list; 1 indicates if a target was observed
end

Base.copy(state::State3d) = State3d(deepcopy(state.koe), deepcopy(state.attitude), deepcopy(state.dt), deepcopy(state.target_list), deepcopy(state.observed_list))

function create_target_list(csv_path)
    # all_data = CSV.read(csv_path, DataFrame)
    all_data = CSV.read(csv_path, DataFrame)
    # all_data.observed .= 0  
    target_list = Tuple.(eachrow(all_data))
    return target_list, length(target_list)
end

function create_target_list_3d(csv_path)
    # all_data = CSV.read(csv_path, DataFrame)
    all_data = CSV.read(csv_path, DataFrame)
    # all_data.observed .= 0  
    target_list = Tuple.(eachrow(all_data[:,3:end]))
    return target_list, length(target_list)
end

""" 
Action space:
1 -> Do nothing OR: rotate to neutral/some intermediate position
2 -> Image target 1
.
.
.
n+1 -> Image target n 
"""

function TR(s, a)

    max_x_ang = 15 
    max_y_ang = 15

    y_noise_mag = .01 
    x_noise_mag = .01

    x_dist = Normal(0, x_noise_mag)
    y_dist = Normal(0, y_noise_mag)

    obs_list = deepcopy(s.observed_list)

    if a == 1
        # No changes to rewards, etc
        # println("Action 1: do nothing")

        R = 0
        attitude = s.attitude

    else
        target = s.target_list[a-1]
        # calculate the required slew angle
        x_angle = atand((target[1] - s.x) / s.alt) # not worried about using atan2 since these are smallish angles
        # println(a[1])
        y_angle = atand((target[2] - s.y) / s.alt)
        # println(a[2])

        if (abs(x_angle) > max_x_ang) || (abs(y_angle) > max_y_ang) 
            # exceeded maximum angle

            # println("Out of slew range. Doing nothing.")
            R = 0
            attitude = s.attitude

        elseif s.observed_list[a-1] == 1
            # already observed
            # we may never reach this state, but including for now just in case

            # println("Already observed this target. Doing nothing.")
            R = 0
            attitude = s.attitude

        else
            # this action is allowed

            attitude = (x_angle, y_angle)    # new attitude after imaging target
            R = target[3]    # getting the reward
            obs_list[a-1] = 1 # set this to 1 to flag that it's been observed

            # println("Slewing to target and imaging.")

        end
    end

    # no matter the action, we continue down the orbit track
    # for temporary test case: 
    y_ = s.y + dydt + rand(y_dist)
    x_ = s.x + rand(x_dist)

    # println("Next time step")

    # for final version: this is where we propagate dynamics TODO

    # TODO: update the list of available targets and actions based on visible horizon

    s_new = State(x_,y_,s.dydt,s.alt, attitude, s.target_list,obs_list) 

    return s_new, R
end

function get_slew_angle(koe, target_tup, dt_JD)
    # Technically this is just getting the angle to the target in RTN frame, NOT the slew angle from current attitude
    mu = 3.986004418e5

    # target tuple is ECEF -- need to convert to ECI
    target_pos = ECEF_to_ECI([target_tup[1], target_tup[2], target_tup[3], 0, 0, 0], dt_JD)
    print("Target position in ECI: ")
    println(target_pos)

    observer_pos = koe2cart(koe, mu)
    print("Observer position in ECI: ")
    println(observer_pos)

    R_eci2rtn = ECI_to_RTN_matrix(observer_pos)

    # println(target_pos)
    # println(observer_pos)

    look_vec_rtn = R_eci2rtn * (target_pos[1:3] .- observer_pos[1:3]) 
    print("Look vector in RTN: ")
    println(look_vec_rtn)

    # get the T and N angles by flattening the look vector into their planes
    T_ang = 90 - acosd(look_vec_rtn[2] / norm( look_vec_rtn[1:2] ))

    N_ang = 90 - acosd(look_vec_rtn[3] / norm( [look_vec_rtn[1], look_vec_rtn[3]] ))

    return (N_ang, T_ang) # Cross-track, along-track

end

function ECI_to_RTN_matrix(rv)
    # Get the RV frame (aka RTN or LVLH)
    r_u = (rv[1:3] / norm(rv[1:3])) # unit vector pointing to nadir
    t_u = (rv[4:6] / norm(rv[4:6])) # unit vector pointing in velocity direction, i.e. along track
    n_u = cross(t_u, r_u) # unit vector in cross-track direction, pointing to the left

    # using method adapted from Duncan Eddy's SatelliteDynamics.jl -- need to cite in report
    # https://github.com/sisl/SatelliteDynamics.jl/blob/46f6c9265b1e648dd3891ad593b122a5d0bfa908/src/reference_systems.jl
    # Note: since we assume a perfectly circular orbit, we can simplify and use v direction as T vec

    R_mat = hcat(r_u, t_u, n_u)
    return R_mat
end

function TR_orbit(s, a)
    # if no time step specified, use 1
    return TR_orbit(s, a, 1)
end

function TR_orbit(s, a, time_step)

    max_t_ang = 15  # max slew angle - in/along track
    max_c_ang = 15 # max slew angle - out of track
    mu = 3.986004418e5 #km^3/m^2   
    Re = 6371 # earth radius, km   # TODO confirm this matches what's in other funcs
    fov = 5 # degrees from boresight

    c_error_mag = 1 # defining error in pointing angle after slew, degrees
    t_error_mag = 1

    c_dist = Normal(0, c_error_mag)
    t_dist = Normal(0, t_error_mag)


    obs_list = deepcopy(s.observed_list)
    # println(length(obs_list))

    rv = koe2cart(s.koe, mu)
    print("Observer state: ")
    println(rv)
    r_u = (rv[1:3] / norm(rv[1:3])) # unit vector pointing to nadir
    t_u = (rv[4:6] / norm(rv[4:6])) # unit vector pointing in velocity direction, i.e. along track
    n_u = cross(t_u, r_u) # unit vector in cross-track direction, pointing to the left

    if a == 1
        # No changes to rewards, etc
        # println("Action 1: do nothing")

        R = 0
        attitude = s.attitude

    else
        target = s.target_list[a-1]
        target_dist = Normal(target[4], target[5])
        # calculate the required slew angle
        # convert target pos to ECI
                
        angs = get_slew_angle(s.koe, target, s.dt)
        print("Target angle: ")
        println(angs)

        # calculate the visible horizon angle and distance
        horizon_angle = acosd(Re/norm(rv[1:3]))
        horizon_dist = sqrt( norm(rv[1:3])^2 - Re^2 )
        println(horizon_angle)
        println(horizon_dist)

        # Check for constraint violations
        out_of_range_c = abs(angs[1] - s.attitude[1]) > max_c_ang
        out_of_range_t = abs(angs[2] - s.attitude[2]) > max_t_ang
        beyond_horizon_c = abs(angs[1]) > horizon_angle
        beyond_horizon_t = abs(angs[2]) > horizon_angle

        target_pos = ECEF_to_ECI([target[1], target[2], target[3], 0, 0, 0], s.dt)
        print("Target state: ")
        println(target_pos)
        observer_pos = koe2cart(koe, mu)
        R_eci2rtn = ECI_to_RTN_matrix(observer_pos)
        look_vec_rtn = R_eci2rtn * (target_pos[1:3] .- observer_pos[1:3]) 
        far_side = norm(look_vec_rtn) > horizon_dist
        println(norm(look_vec_rtn))
        println(look_vec_rtn)

        if obs_list[a-1] == 1
            # already observed

            println("Already observed this target. Doing nothing.")
            R = 0 
            attitude = s.attitude
        elseif far_side
            # we definitely can't see it if it's on the far side of the earth
            println("Target is beyond line of sight. Doing nothing.")
            R = 0
            attitude = s.attitude

        elseif out_of_range_t || out_of_range_c || beyond_horizon_t || beyond_horizon_c 
            # exceeded maximum angle

            # we can slew only as far as the maximum slew angle / horizon_angle
            attitude = s.attitude # just initializing as something
            if angs[1] > (attitude[1] + max_c_ang)
                # we have exceeded the range in + direction
                attitude[1] = minimum([ attitude[1] + max_c_ang, horizon_angle ]) + rand(c_dist) # choose whatever is lowest
            else
                # we have exceeded the range in - direction
                attitude[1] = maximum([ attitude[1] - max_c_ang, -horizon_angle ]) + rand(c_dist) # choose whatever is highest (closer to 0)
            end

            if angs[2] > (attitude[2] + max_t_ang)
                # we have exceeded the range in + direction
                attitude[2] = minimum([ attitude[2] + max_t_ang, horizon_angle ]) + rand(t_dist) # choose whatever is lowest
            else
                # we have exceeded the range in - direction
                attitude[2] = maximum([ attitude[2] - max_t_ang, -horizon_angle ]) + rand(t_dist) # choose whatever is highest (closer tp 0)
            end

            println("Out of slew range. Slewed to limit of slew range (with error).")

            if (angs[1] >= attitude[1] - fov) & (angs[1] <= attitude[1] + fov) &
                 (angs[2] >= attitude[2] - fov) & (angs[2] <= attitude[2] + fov)
                # target is in the field of view
                println("Imaged target anyway!")
                # println(target_dist)
                # println(rand(target_dist))
                R = maximum( [rand(target_dist), 0 ])
                obs_list[a-1] = 1 # set this to 1 to flag that it's been observed

            else
                # target outside field of view
                println("Did not successfully image target.")
                R = 0

            end

        else
            # this action is within constraints --> allowed


            attitude = state.attitude # just to initialize it
            attitude[1] = angs[1] + rand(c_dist) 
            attitude[2] = angs[2] + rand(t_dist) 

            println("Slewing to target.")

            if (angs[1] >= attitude[1] - fov) & (angs[1] <= attitude[1] + fov) &
                 (angs[2] >= attitude[2] - fov) & (angs[2] <= attitude[2] + fov)
                # target is in the field of view
                println("Target in FOV - imaged!")
                R = maximum( [rand(target_dist), 0] )
                obs_list[a-1] = 1 # set this to 1 to flag that it's been observed

            else
                # target outside field of view
                println("Target out of FOV - did not image.")
                R = 0

            end


        end
    end

    # no matter the action, we continue down the orbit track
    new_koe = kepler_dyn(copy(s.koe), time_step, mu )
    dt = s.dt + time_step

    println("Next time step")

    s_new = State3d(new_koe, attitude, dt, s.target_list,obs_list) 
 
    return s_new, R
end