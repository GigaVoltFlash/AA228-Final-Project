
# koe = [a,e,i,Ω,ω,M]
function kepler_dyn(koe::Vector, dt, param)

    a = koe[1]
    n = sqrt(param.mu/a^3)
    koe[6] = koe[6] + n*dt

    return koe
end

function cart_dyn(cart::Vector, dt, param)

    koe = cart2koe(cart, param.mu)
    koe = kepler_dyn(koe, dt, param)
    cart = koe2cart(koe, param.mu)

    return cart
end


function cart2koe(cart, mu)
    r = norm(cart[1:3])
    v = norm(cart[4:6])

    ene = v^2/2 - mu/r
    a = -mu/(2*ene)
    hvec = cross(cart[1:3], cart[4:6])
    evec = 1/mu * cross(cart[4:6], hvec) - cart[1:3]/r
    zvec = [0,0,1]
    nvec = cross(zvec, hvec) / norm(cross(zvec, hvec))  # line of node

    i = acos(dot(zvec, hvec)/norm(hvec))  # 0<i<π
    Ω = atan2(nvec[2], nvec[1])  # 0<Ω<2π (no ambiguity)
    ω = acos(dot(nvec, evec)/norm(evec))  
    if evec[3] < 0
        ω = 2*pi - ω
    end
    f = acos(dot(evec, cart[1:3])/r/norm(evec))  
    if dot(cart[1:3], cart[4:6]) < 0
        f = 2*pi - f
    end

    E = atan2(sqrt(1-e^2)*sin(f), e+cos(f))   # no ambiguity
    M = E - e*sin(E)  # no ambiguity

    return [a,e,i,Ω,ω,M]
end 


# koe = [a,e,i,Ω,ω,M]
function koe2cart(koe, mu)
    a = koe[1]
    e = koe[2]
    i = koe[3]
    Ω = koe[4]
    ω = koe[5]
    M = koe[6]
    
    E = solve_kepler(M, e)
    f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2))

    r = a*(1-e^2)/(1+e*cos(f))
    p = a*(1-e^2)

    # perifocal frame
    rvec = [r*cos(f), r*sin(f), 0]  
    vvec = sqrt(mu/p)*[-sin(f), e+cos(f), 0]

    # rotation matrices
    r = Rz(ω) * Rx(i) * Rz(Ω) * rvec
    v = Rz(ω) * Rx(i) * Rz(Ω) * vvec 

    return cat(r, v)...  # concatenate r and v
end


function solve_kepler(M, e)
    E = M
    for _ = 1:100
        E_next = E - (E - e*sin(E) - M)/(1 - e*cos(E))

        if abs(E_next - E) < 1e-5
            break
        end
    end

    return E
end


# rotation matrix: passive rotation for coordinate transformation 
function Rz(x)
    return [cos(x) sin(x) 0; -sin(x) cos(x) 0; 0 0 1]
end 

function Rx(x)
    return [1 0 0; 0 cos(x) sin(x); 0 -sin(x) cos(x)]
end




function TR(s, a)

    newstate = s 

    max_x_ang = 30 
    max_y_ang = 30

    y_noise_mag = .01 
    x_noise_mag = .01

    x_dist = Normal(0, x_noise_mag)
    y_dist = Normal(0, y_noise_mag)

    obs_list = np.zeros(size(s.obs_list,1))

    if a == 1
        # No changes to rewards, etc
        # println("Action 1: do nothing")

        R = 0

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

        elseif s.observed_list[a-1] == 1
            # already observed
            # we may never reach this state, but including for now just in case

            # println("Already observed this target. Doing nothing.")
            R = 0

        else
            # this action is allowed

            attitude = (x_angle, y_angle)    # new attitude after imaging target
            R = target[3]    # getting the reward
            # state.target_list[a-1][3] = 0   # set the reward in the image tuple to 0 so we don't image it again
            terget_list[a-1] = 1 # set this to 1 to flag that it's been observed

            # println("Slewing to target and imaging.")

        end
    end

    # no matter the action, we continue down the orbit track
    # for temporary test case: 
    y_ += dydt + rand(y_dist)
    x_ += rand(x_dist)

    # println("Next time step")

    # for final version: this is where we propagate dynamics TODO

    # TODO: update the list of available targets and actions based on visible horizon

    s_new = State(x,y,s.dydt,s.alt,s.target_list,obs_list) 

    return s_new, R
end

