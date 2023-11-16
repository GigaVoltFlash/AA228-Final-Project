
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