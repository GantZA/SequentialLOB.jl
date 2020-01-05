mutable struct ObjectiveSurface
    param1_name::String
    param1_values::Array{Any,1}
    param2_name::String
    param2_values::Array{Any,1}
    surface_points::Int64
    replications::Int64
    params::Dict{String, Any}
end


function (os::ObjectiveSurface)(seed, para=false)
    iterations = os.surface_points^2

    seeds = Int.(rand(MersenneTwister(seed), UInt32, iterations))
    if para==true
        sample_price_paths = SharedArray{Float64,2}((os.params["T"]+1, iterations*os.replications))
        @sync @distributed for i in 1:iterations
            os.params[os.param1_name] = os.param1_values[i]
            os.params[os.param2_name] = os.param2_values[i]
            rdpp = ReactionDiffusionPricePaths(os.params)
            index_start = (i-1) * os.replications + 1
            index_end = i * os.replications
            _, _, sample_price_paths[:, index_start:index_end], _, _ = rdpp(seeds[i])
        end
        return sample_price_paths
    else
        sample_price_paths = Array{Float64,2}(undef, os.params["T"]+1, iterations*os.replications)
        for i in 1:iterations
            os.params[os.param1_name] = os.param1_values[i]
            os.params[os.param2_name] = os.param2_values[i]
            rdpp = ReactionDiffusionPricePaths(os.params)
            index_start = (i-1) * os.replications + 1
            index_end = i * os.replications
            _, _, sample_price_paths[:, index_start:index_end], _, _ = rdpp(seeds[i])
        end
        return sample_price_paths
    end
end
