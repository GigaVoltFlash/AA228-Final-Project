include("state.jl")
if false

    # Create an intial state
    x0 = 0 # checl: are the test points centered around x=0?
    y0 = 0
    dydt = 1
    alt = 50
    attitude = [0.0, 0.0]
    target_list, n_targets = create_target_list("src/obs_site.csv")

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
    koe = [6571, 0.001, 0, 0, 0, 0]
    att = [0,0]
    dt = 0
    target_list = [(6371, 0, 0, 100), (0, 6371, 0, 100)]
    n_targets = 2
    observed_list = [0, 0]
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

    reward_total = 0
    for t in 0:10
        global state, reward_total
        # a = rand(A)
        a = 2
        print("Choosing action ")
        print(a)
        print("\n")

        state, reward = TR_orbit(state, a)

        reward_total += reward

    end

    print("Total reward: ")
    println(reward_total)







end