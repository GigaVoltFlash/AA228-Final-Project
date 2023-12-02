using CSV, Random, DataFrames
using GLMakie
using GeoDatasets

function create_obs_sites_csv(x_max, y_max, num_rows::Int)
    # Create a DataFrame with specified columns
    df = DataFrame(x = round.(-x_max .+ (2*x_max) .* rand(num_rows), digits=2), 
                   y = round.(y_max .* rand(num_rows), digits=2), 
                   r = round.(rand(num_rows), digits=2))
    
    # Save the DataFrame to a CSV file
end


function map_on_earth(λ0, i, Δλ, Δφ, num_site::Int)
    """
        inputs: 
            λ0: longitude of the center of the target map (equatorial)
            i: inclination of the orbit
            Δλ: longitude span of the map
            Δφ: latitude span of the map
            num_site: number of observation sites
    """    
    R_E = 6375  # km

    df = DataFrame(λ = round.(λ0 .+ Δλ .* (rand(num_site).*2 .- 1), digits=2), 
                   φ = round.(Δφ .* (rand(num_site).*2 .- 1), digits=2))
    df.λ = df.λ - df.φ .* tand(i-90)

    df.x = R_E .* cosd.(df.φ) .* cosd.(df.λ)
    df.y = R_E .* cosd.(df.φ) .* sind.(df.λ)
    df.z = R_E .* sind.(df.φ)

    return df 
end

lon,lat,data = GeoDatasets.landseamask(;resolution='c',grid=5)

df = map_on_earth(-110, 98, 10, 60, 100)
df = CSV.write("obs_site_Earth.csv", df)
f = Figure(size = (1200, 800))
ax = Axis(f[1, 1], aspect=DataAspect(), title = "Target sites")
scatter!(ax, df.λ, df.φ, overdraw = true, color = :red, markersize = 3)
contour!(ax, lon, lat, data, levels=[0.5], color=:black, linewidth=0.5)
f

