using DataFrames
using CSV

mutable struct State
    # lat             # in degrees
    # lon             # in degrees
    x               # for the simple test case only - in km
    y               # for the simple test case only - in km
    dydt            # "velocity" of the s/c -- how far does it move in a time step. Probably can be defined outside this struct
    alt             # in km(?)
    attitude        # Tuple of angles - for now: (x_angle, y_angle)
    target_list     # n entries in list; each entry is a tuple of (x, y, r)
end

function create_target_list(csv_path)
    # all_data = CSV.read(csv_path, DataFrame)
    all_data = CSV.read(csv_path, DataFrame)
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

function transition(state, a)

    max_x_ang = 30 
    max_y_ang = 60

    if a == 1
        # No changes to rewards, etc
        println("Action 1: do nothing")

        # TODO figure out if we want to do anything here

    else
        target = state.target_list[a-1]
        # calculate the required slew angle
        x_angle = atand((target[1] - state.x) / state.alt) # not worried about using atan2 since these are smallish angles
        # println(a[1])
        y_angle = atand((target[2] - state.y) / state.alt)
        # println(a[2])

        if (abs(x_angle) > max_x_ang) || (abs(y_angle) > max_y_ang) 
            # exceeded maximum angle

            # TODO: figure out what happens here
            # for now: do nothing

            println("Out of slew range. Doing nothing.")

        else
            # this action is allowed

            # TODO collect reward
            
            # TODO Somewhere, we need to remove the action from the list (or otherwise mark that it was completed) if it's actually taken

            state.attitude = (x_angle, y_angle)

            println("Slewing to target and imaging.")

        end

    end

    # no matter the action, we continue down the orbit track
    # for temporary test case: 
    state.y += dydt

    println("Next time step")

    # for final version: this is where we propagate dynamics TODO

end