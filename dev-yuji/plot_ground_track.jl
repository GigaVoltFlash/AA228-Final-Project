using CSV, Random, DataFrames
using GLMakie
using GeoDatasets
include("../src/orbits.jl")

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

mu_E = 398600.4418  # [km^3/s^2]


koe0 = [
    7057.0,  # a [km]
    0.000879,  # e 
    deg2rad(98.12),  # i
    deg2rad(255.09),  # Ω [rad]
    deg2rad(225.37),  # ω [rad]
    deg2rad(75.0),  # M [rad]
]

# oe propagation
dt = 20 
time = 0:20:2000
state = zeros(6, length(time))
state_ECEF = zeros(6, length(time))
state_ECI = zeros(6, length(time))
lat_long = zeros(2, length(time))

for (i, t) in enumerate(time)
    koe = kepler_dyn(koe0, dt, mu_E)
    rv = koe2cart(koe, mu_E)
    state[:,i] = rv
    state_ECEF[:,i] = ECI_to_ECEF(rv, t)
    λ = atand(state_ECEF[2,i], state_ECEF[1,i])
    φ = atand(state_ECEF[3,i], sqrt(state_ECEF[1,i]^2 + state_ECEF[2,i]^2))
    lat_long[:,i]  = [λ; φ]
end

# load df 
df = CSV.read("src/obs_site_Earth2.csv", DataFrame)

# plotting 
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
f = Figure(size = (1200, 800))
ax = Axis(f[1, 1], aspect=DataAspect(), title = "Target sites")
scatter!(ax, df.lambda, df.phi, overdraw = true, color = :red, markersize = 3)
contour!(ax, lon, lat, data, levels=[0.5], color=:black, linewidth=0.5)
lines!(ax, lat_long[1,:], lat_long[2,:], overdraw = true, color = :blue, linewidth = 2,)
ylims!(ax, low = -90, high = 90)
xlims!(ax, low = -180, high = 180)
savefig(f, "random_targeting.png")
