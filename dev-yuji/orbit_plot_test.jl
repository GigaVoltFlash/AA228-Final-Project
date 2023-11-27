using GLMakie 
include("orbits.jl")

mu_E = 398600.4418  # [km^3/s^2]
r_E = 6375.0  # [km]

# Molniya orbit
koe0 = [
    26600.0,  # a [km]
    0.74,  # e 
    deg2rad(63.4),  # i
    0.0,  # Ω [rad]
    deg2rad(270),  # ω [rad]
    0.0,  # M [rad]
]

dt = 20 
time = 0:20:86400 / 4
state = zeros(6, length(time))

for (i, t) in enumerate(time)
    koe = kepler_dyn(koe0, dt, mu_E)
    rv = koe2cart(koe, mu_E)
    state[:,i] = rv
end


# x, y = collect(-8:0.5:8), collect(-8:0.5:8)
# z = [sinc(√(X^2 + Y^2) / π) for X ∈ x, Y ∈ y]
# wireframe(x, y, z, axis=(type=Axis3,), color=:black)

n = 33
u = range(0,stop=2*π,length=n);
v = range(0,stop=π,length=n);
x = zeros(n,n); y = zeros(n,n); z = zeros(n,n)
for i in 1:n
    for j in 1:n
        x[i,j] = r_E*cos.(u[i]) * sin(v[j]);
        y[i,j] = r_E*sin.(u[i]) * sin(v[j]);
        z[i,j] = r_E*cos(v[j]);
    end
end

f = Figure()
ax = Axis3(f[1, 1], aspect = :data,     title = "Molniya orbit (ECI)",)
# println("orbit: ",  state[1,:], state[2,:], state[3,:])
wireframe!(ax, x,y,z, overdraw = true, transparency = true, color = (:black, 0.1))
lines!(ax, state[1,:], state[2,:], state[3,:], overdraw = true, color = :red, linewidth = 2, )
# scatter!(ax, state[1,:], state[2,:], state[3,:], overdraw = true, color = :red, markersize = 3)
f