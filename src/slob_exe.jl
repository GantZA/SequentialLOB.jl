module SLOBExec

using SequentialLOB

Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    parsed_args = parse_commandline()
    slob_model = SLOB(
        1,
        parsed_args["T"],
        parsed_args["p₀"],
        parsed_args["M"],
        parsed_args["L"],
        parsed_args["D"],
        parsed_args["σ"],
        parsed_args["nu"],
        parsed_args["α"],
        SourceTerm(parsed_args["λ"], parsed_args["μ"]),
    )

    mid_price_paths = slob_model(parsed_args["SEED"])
    print(mid_price_paths[:, 1])
    return 0
end


end
