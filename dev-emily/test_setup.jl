include("state.jl")
if false

    # Create an intial state
    x0 = 0 # checl: are the test points centered around x=0?
    y0 = 0
    dydt = 1
    alt = 50
    attitude = [0.0, 0.0]
    target_list, n_targets = create_target_list_3d("src/obs_site.csv")

    state = State(x0, y0, dydt, alt, attitude, target_list, zeros(n_targets))
    # print(state)
    # create action list
    """ 
    Action space:
    1 -> Do nothing OR: rotate to neutral/some intermediate position
    2 -> Image target 1
    .
    .
    .
    n+1 -> Image target n 
    """
    A = 1:n_targets + 1

    for t = 0:60
        global state
        a = rand(A)
        print("Choosing action ")
        print(a)
        print("\n")
        state, reward = TR(state, a)
        # println(state)

    end
else
    # create initial 3d state
    # koe = [6571, 0.001, 0, 0, 0, 0]
    
    koe = [
    7057.0,  # a [km]
    0.000879,  # e 
    deg2rad(98.12),  # i
    deg2rad(260.09),  # Ω [rad]
    deg2rad(221.37),  # ω [rad]
    0.0,  # M [rad]
    ]
    att = [0.0,0.0]
    dt = 0
    target_list, n_targets = create_target_list_3d("src/obs_site_Earth.csv")
    observed_list = zeros(n_targets)
    state = State3d(koe, att, dt, target_list, observed_list)



    """ 
    Action space:
    1 -> Do nothing OR: rotate to neutral/some intermediate position
    2 -> Image target 1
    .
    .
    .
    n+1 -> Image target n 
    """
    A = 1:n_targets + 1

    time_step = 30 # seconds
    reward_total = 0
    for t in 0:time_step:30*60
        global state, reward_total
        a = rand(A)
        # a = 92
        print("Choosing action ")
        print(a)
        print("\n")

        state, reward = TR_orbit(state, a, time_step)

        reward_total += reward

    end

    print("Total reward: ")
    println(reward_total)



end