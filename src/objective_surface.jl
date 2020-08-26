mutable struct ObjectiveSurface
    param1_name::String
    param1_values::Array{Float64,1}
    param2_name::String
    param2_values::Array{Float64,1}
    surface_points::Int64
    replications::Int64
    params::Dict{String, Any}
    weight_matrix::Any
end


function objective_function(parameter_dict::Dict, os::ObjectiveSurface, seed, ind)
    try
        parameter_dict[os.param1_name] = os.param1_values[ind]
        parameter_dict[os.param2_name] = os.param2_values[ind]
        log_price_paths = log.(SLOB(parameter_dict)(seed))
        log_price_paths = vcat(log_price_paths[1:1,:], log_price_paths)
        log_returns = diff(log_price_paths, dims=1)
        objective_value = weighted_moment_distance(os.weight_matrix,
            log_returns)
        return objective_value
    catch e
        return NaN
    end
end


function (os::ObjectiveSurface)(seed, parallel=false)
    iterations = os.surface_points^2

    seeds = Int.(rand(MersenneTwister(seed), UInt32, iterations))
    parameter_vector = fill(os.params, iterations)
    if parallel==true
        os_values = SharedArray{Float64,1}(iterations)
        @sync @distributed for i in 1:iterations
            os_values[i] = objective_function(parameter_vector[i], os, seeds[i], i)
        end
        return Array(os_values)
    else
        os_values = Array{Float64,1}(undef, iterations)
        @showprogress 1 "Iteration: " for i in 1:iterations
            os_values[i] = objective_function(parameter_vector[i], os, seeds[i], i)
        end
        return os_values
    end
end
