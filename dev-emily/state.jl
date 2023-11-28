using DataFrames
using CSV
using Distributions
using LinearAlgebra

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
    attitude        # Tuple of angles - for now: (E_angle, N_angle)
    target_list     # n entries in list; each entry is a tuple of (x, y, z, r) w/ posiitons in ECEF
    observed_list   # n entries in list; 1 indicates if a target was observed
end

Base.copy(state::State3d) = State3d(deepcopy(state.koe), deppcopy(state.attitude), deepcopy(state.target_list), deepcopy(state.observed_list))

function create_target_list(csv_path)
    # all_data = CSV.read(csv_path, DataFrame)
    all_data = CSV.read(csv_path, DataFrame)
    # all_data.observed .= 0  
    target_list = Tuple.(eachrow(all_data))
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

function get_slew_angle(koe, target_tup)
    # target tuple is ECEF






end


function TR_orbit(s, a)

    max_it_ang = 15  # in/along track
    max_oot_ang = 15 # out of track
    mu = 3.986004418e5

    # x_noise_mag = .01
    # y_noise_mag = .01

    # x_dist = Normal(0, x_noise_mag)
    # y_dist = Normal(0, y_noise_mag)

    obs_list = deepcopy(s.observed_list)

    rv = koe2cart(s.koe, mu)
    nadir_u = -(rv[1,:] / norm(rv[1,:])) # unit vector pointing to nadir
    v_u = (rv[2,:] / norm(rv[2,:])) # unit vector pointing in velocity direction, i.e. along track
    ct_u = - cross(v_u, nadir_u) # unit vector in cross-track direction, pointing to the right


    if a == 1
        # No changes to rewards, etc
        # println("Action 1: do nothing")

        R = 0
        attitude = s.attitude

    else
        target = s.target_list[a-1]
        # calculate the required slew angle

        # TODO: Convert this to use ECEF XYZ positions for both
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