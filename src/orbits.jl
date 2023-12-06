using Dates
using LinearAlgebra

# constants 
Ωdot = 2*pi / (23*3600 + 56*60 + 4.09053) 
tjd0 = DateTime(2000,1,1,12,0,0)  # J2000 epoch

"""
Nominal Orbit of the satellite (GeoEye-1 / MAXER)
a = 7,057 km (4,385 mi)
e = 0.000879
i = 98.12 deg. 
Ω = 347.09 deg.  (J2 will rotate so we can choose whatever we want)
ω = 221.37 deg.
Period	98.34 minutes
RAAN	347.09 degrees
"""
        

# koe = [a,e,i,Ω,ω,M]
function kepler_dyn(koe::Vector, dt, mu_E)

    a = koe[1]
    n = sqrt(mu_E/a^3)
    koe_ = copy(koe)
    koe_[6] = koe[6] + n*dt

    return koe_
end

function kepler_dyn!(koe::Vector, dt, mu_E)
    a = koe[1]
    n = sqrt(mu_E/a^3)
    koe[6] = koe[6] + n*dt
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
    r = Rz(-Ω) * Rx(-i) * Rz(-ω) * rvec
    v = Rz(-Ω) * Rx(-i) * Rz(-ω) * vvec 

    return vcat(r, v)  # concatenate r and v
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


function ECI_to_ECEF(rv, dt_JD)
    Ω = Ωdot * dt_JD;
    r = Rz(Ω) *  rv[1:3] 
    v = Rz(Ω) * rv[4:6] - cross([0,0,Ωdot], rv[1:3])
    return vcat(r, v)
end 

function ECEF_to_ECI(rv, dt_JD)
    Ω = Ωdot * dt_JD;
    r = Rz(-Ω) * rv[1:3] 
    v = Rz(-Ω) * rv[4:6] + cross([0,0,Ωdot], rv[1:3])
    return vcat(r, v)
end

# rotation matrix: passive rotation for coordinate transformation 
function Rx(x)
    return [1 0 0; 0 cos(x) sin(x); 0 -sin(x) cos(x)]
end

function Ry(x)
    return [cos(x) 0 -sin(x); 0 1 0; sin(x) 0 cos(x)]
end

function Rz(x)
    return [cos(x) sin(x) 0; -sin(x) cos(x) 0; 0 0 1]
end 



