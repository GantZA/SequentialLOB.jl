# Generate a time series of stock prices by repeatedly solving the SPDE
# from Mastromatteo et al. (2014) to model the latent order book.
# using the previously generated mid price as the new initial mid price and
# boundary mid point during each iteration to generate a sequence of T prices,
# where T is the number of measured prices stored in the loaded
# ObjectiveFunction object.

mutable struct ReactionDiffusionPricePaths
    num_paths::Int64
    T::Int64
    p₀::Float64
    M::Int64
    β::Float64
    L::Float64
    D::Float64
    σ::Float64
    nu::Float64
    α::Float64
    source_term::SourceTerm
    x::Array{Float64, 1}
    Δx::Float64
    Δt::Float64
end
function ReactionDiffusionPricePaths(num_paths::Int64, T::Int64, p₀::Float64,
    M::Int64, β::Float64, L::Real, D::Float64, σ::Float64, nu::Float64,
    α::Float64, source_term::SourceTerm)
    logger = FileLogger(Dict(Logging.Info => "info.log", Logging.Error => "error.log"), append=false)
    global oldglobal = global_logger(logger)
    x₀ = p₀ - 0.5*L
    xₘ = p₀ + 0.5*L
    @assert x₀ >= 0
    # if D < 1.5 * σ
    #     @warn "D=$D is less than 1.5 times σ=$σ which runs the risk of unstable LOB ghost points"
    # end
    # Δx = L/M
    # if (β*Δx/D < log(1/3)/2) | (β*Δx/D > log(3)/2)
    #     @warn "The choice of β, D, L, M might result in jump probabilities > 0.75, thus causing large shifts in the LOB density"
    # end

    x = collect(Float64, range(x₀, xₘ, length=M+1))
    Δx = L/M
    Δt = (Δx^2) / (2.0*D)
    return ReactionDiffusionPricePaths(num_paths, T, p₀, M, β,
        L, D, σ, nu, α, source_term, x, Δx, Δt)
end


ReactionDiffusionPricePaths(dict)=ReactionDiffusionPricePaths(
    dict["num_paths"], dict["T"], dict["p₀"], dict["M"],
    dict["β"], dict["L"], dict["D"],
    dict["σ"], dict["nu"], dict["α"], SourceTerm(dict["λ"], dict["μ"]))


ReactionDiffusionPricePaths(;num_paths::Int64=1, T::Int64=100, p₀::Real=100.0,
    M::Int64=100, β::Real=1.0, L::Real=50.0, D::Real=4.0, σ::Real=1.0,
    nu::Real=0.0, α::Real=1.0, λ::Real=1.0, μ::Real=0.5) =
    ReactionDiffusionPricePaths(num_paths, T, p₀, M, β, L, D, σ, nu, α,
    SourceTerm(λ, μ))


function (rdpp::ReactionDiffusionPricePaths)(seed::Int=-1)
    if seed == -1
            seeds = Int.(rand(MersenneTwister(), UInt32, rdpp.num_paths))
        else
            seeds = Int.(rand(MersenneTwister(seed), UInt32, rdpp.num_paths))
    end
    time_steps = floor(Int, rdpp.T / rdpp.Δt)

    raw_price_paths = ones(Float64, time_steps + 1, rdpp.num_paths)
    raw_price_paths[1, :] .= rdpp.p₀

    sample_price_paths = ones(Float64, rdpp.T + 1, rdpp.num_paths)
    sample_price_paths[1, :] .= rdpp.p₀

    lob_densities = zeros(Float64, rdpp.M+1, time_steps + 1, rdpp.num_paths)
    P⁺s = ones(Float64, time_steps, rdpp.num_paths)
    P⁻s = ones(Float64, time_steps, rdpp.num_paths)
    Ps = ones(Float64, time_steps, rdpp.num_paths)
    for path in 1:rdpp.num_paths
        Random.seed!(seeds[path])
        @info "path $path with seed $(seeds[path])"
        lob_densities[:, :, path], raw_price_paths[:, path],
            sample_price_paths[:, path], P⁺s[:, path],
            P⁻s[:, path], Ps[:, path]  = dtrw_solver(rdpp)
    end

    return lob_densities, raw_price_paths, sample_price_paths, P⁺s, P⁻s, Ps
end
