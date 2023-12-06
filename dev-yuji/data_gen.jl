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

    df = DataFrame(lambda = round.(λ0 .+ Δλ .* (rand(num_site).*2 .- 1), digits=2), 
                   phi = round.(Δφ .* (rand(num_site).*2 .- 1), digits=2))
    df.lambda = df.lambda - df.phi .* tand(i-90)

    df.x = R_E .* cosd.(df.phi) .* cosd.(df.lambda)
    df.y = R_E .* cosd.(df.phi) .* sind.(df.lambda)
    df.z = R_E .* sind.(df.phi)
    df.r_mean = rand(num_site)
    df.r_std  = rand(num_site) * 0.3

    return df 
end


df = map_on_earth(-110, 100, 10, 60, 50)
df = CSV.write("obs_site_Earth_50.csv", df)


