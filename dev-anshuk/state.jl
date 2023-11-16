using DataFrames
using CSV
using Distributions

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