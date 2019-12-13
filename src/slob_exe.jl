module LOBMExec

using LatentOrderBookModel

Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    parsed_args = parse_commandline()
    rdpp = ReactionDiffusionPricePaths(
        parsed_args["num_paths"],
        parsed_args["T"],
        parsed_args["p₀"],
        parsed_args["M"],
        parsed_args["β"],
        parsed_args["L"],
        parsed_args["D"],
        parsed_args["σ"],
        parsed_args["nu"],
        parsed_args["α"],
        SourceTerm(parsed_args["λ"], parsed_args["μ"]),
    )

    lob_densities,
    price_paths,
    mid_price_bars,
    P⁺s,
    P⁻s = rdpp(parsed_args["SEED"])
    print(price_paths[:, 1])
    return 0
end


end
