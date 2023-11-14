using CSV, Random, DataFrames

function create_obs_sites_csv(fname, x_max, y_max, num_rows::Int)
    # Create a DataFrame with specified columns
    df = DataFrame(x = round.(-x_max .+ (2*x_max) .* rand(num_rows), digits=2), 
                   y = round.(-y_max .+ (2*y_max) .* rand(num_rows), digits=2), 
                   r = round.(rand(num_rows), digits=2))
    
    # Save the DataFrame to a CSV file
    CSV.write(fname, df)
end

create_obs_sites_csv("obs_site.csv", 20, 60, 100)