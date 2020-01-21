mutable struct ObjectiveSurface
    param1_name::String
    param1_values::Array{Float64,1}
    param2_name::String
    param2_values::Array{Float64,1}
    surface_points::Int64
    replications::Int64
    params::Dict{String, Any}
    weight_matrix::BlockBootstrapWeightMatrix
end


function objective_function(os::ObjectiveSurface, seed)
    try
        slob = SLOB(os.params)
        mid_price_paths = slob(seed)
        objective_value = weighted_moment_distance(os.weight_matrix,
            mid_price_paths)
        return objective_value
    catch e
        return NaN
    end
end


function (os::ObjectiveSurface)(seed, parallel=false)
    iterations = os.surface_points^2

    seeds = Int.(rand(MersenneTwister(seed), UInt32, iterations))
    if parallel==true
        os_values = SharedArray{Float64,1}(iterations*os.replications)
        @sync @distributed for i in 1:iterations
            os.params[os.param1_name] = os.param1_values[i]
            os.params[os.param2_name] = os.param2_values[i]
            os_values[i] = objective_function(os, seeds[i])
        end
        return os_values
    else
        os_values = Array{Float64,1}(undef, iterations*os.replications)
        for i in 1:iterations
            os.params[os.param1_name] = os.param1_values[i]
            os.params[os.param2_name] = os.param2_values[i]
            os_values[i] = objective_function(os, seeds[i])
        end
        return os_values
    end
end
