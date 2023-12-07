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

koe0_ = copy(koe0)

# oe propagation
dt = 30
time = 0:dt:2000
state = zeros(6, length(time))
state_ECEF = zeros(6, length(time))
state_ECI = zeros(6, length(time))
lat_long = zeros(2, length(time))

for (i, t) in enumerate(time)
    kepler_dyn!(koe0, dt, mu_E)
    rv = koe2cart(koe0, mu_E)
    state[:,i] = rv
    state_ECEF[:,i] = ECI_to_ECEF(rv, t)
    λ = atand(state_ECEF[2,i], state_ECEF[1,i])
    φ = atand(state_ECEF[3,i], sqrt(state_ECEF[1,i]^2 + state_ECEF[2,i]^2))
    lat_long[:,i]  = [λ; φ]
end

# df_traj = (
#     t = time,
#     x_ECI = state[1,:], y_ECI = state[2,:], z_ECI = state[3,:], 
#     vx_ECI = state[4,:], vy_ECI = state[5,:], vz_ECI = state[6,:],
#     x_ECEF = state_ECEF[1,:], y_ECEF = state_ECEF[2,:], z_ECEF = state_ECEF[3,:],
#     vx_ECEF = state_ECEF[4,:], vy_ECEF = state_ECEF[5,:], vz_ECEF = state_ECEF[6,:],
#     lambda = lat_long[1,:], phi = lat_long[2,:] 
#     )
# CSV.write("traj_info.csv", df_traj)

# load df 
df = CSV.read("src/obs_site_Earth2.csv", DataFrame)

# plotting 
fontsize_theme = Theme(fontsize = 20)
set_theme!(fontsize_theme)
lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)
f = Figure(size = (3000, 2400))
ax = Axis(f[1, 1], aspect=DataAspect(), title = " ", titlesize = 20, 
            xlabel = "Longitude [deg.]", ylabel = "Latitude [deg.]", xlabelsize = 20, ylabelsize = 20,
            )
# scatter!(ax, df.lambda, df.phi, overdraw = true, color = :red, markersize = 3)
contour!(ax, lon, lat, data, levels=[0.5], color=:black, linewidth=0.5)
traj1 = lines!(ax, lat_long[1,:], lat_long[2,:], overdraw = true, color = :blue, linewidth = 2,)

k = 10.580 
Δλ = 10
p_br = [-110 + k + Δλ, -60]
p_bl = [-110 + k - Δλ, -60]
p_tr = [-110 - k + Δλ,  60]
p_tl = [-110 - k - Δλ,  60]

l1 = lines!(ax, [p_br[1], p_bl[1]], [p_br[2], p_bl[2]], overdraw=true, color=:orange, linestyle=:solid, linewidth=2,)
l2 = lines!(ax, [p_tr[1], p_tl[1]], [p_tr[2], p_tl[2]], overdraw=true, color=:orange, linestyle=:solid, linewidth=2,)
l3 = lines!(ax, [p_br[1], p_tr[1]], [p_br[2], p_tr[2]], overdraw=true, color=:orange, linestyle=:solid, linewidth=2,)
l4 = lines!(ax, [p_bl[1], p_tl[1]], [p_bl[2], p_tl[2]], overdraw=true, color=:orange, linestyle=:solid, linewidth=2,)


x0 = scatter!(ax, [lat_long[1,1]],   [lat_long[2,1]], overdraw = true, marker=:circle , color = :red, markersize = 10)
xf = scatter!(ax, [lat_long[1,end]], [lat_long[2,end]], overdraw = true, marker=:star5 , color = :purple, markersize = 20)
ylims!(ax, low = -90, high = 90)
xlims!(ax, low = -180, high = 180)
ax.xticks = ([-180, -90, 0, 90, 180])
# savefig(f, "random_targeting.png")

# oe propagation
time = 2000:dt:6000
state = zeros(6, length(time))
state_ECEF = zeros(6, length(time))
state_ECI = zeros(6, length(time))
lat_long = zeros(2, length(time))

for (i, t) in enumerate(time)
    kepler_dyn!(koe0, dt, mu_E)
    rv = koe2cart(koe0, mu_E)
    state[:,i] = rv
    state_ECEF[:,i] = ECI_to_ECEF(rv, t)
    λ = atand(state_ECEF[2,i], state_ECEF[1,i])
    φ = atand(state_ECEF[3,i], sqrt(state_ECEF[1,i]^2 + state_ECEF[2,i]^2))
    lat_long[:,i]  = [λ; φ]
end

k = 
lat_long1 = lat_long[:, 1:12]
lat_long2 = lat_long[:, 13:end]
traj2 = lines!(ax, lat_long1[1,:], lat_long1[2,:], overdraw=true, color=:blue, linestyle=:dash, linewidth=2,)
traj2 = lines!(ax, lat_long2[1,:], lat_long2[2,:], overdraw=true, color=:blue, linestyle=:dash, linewidth=2,)

axislegend(ax, 
    [x0, xf, traj1, traj2, [l1, l2, l3, l4]],
    ["initial state", "terminal state", "ground track (2000s)", "ground track (100 min)", "obs. region", ], 
    position = :rb, labelsize = 15)

save("ground_track.png", f)
f