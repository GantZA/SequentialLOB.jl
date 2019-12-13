using LatentOrderBookModel

function main(output = "stdout")
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
    if output == "stdout"
        lob_densities,
        raw_price_paths,
        sample_price_paths,
        P⁺s,
        P⁻s = rdpp(parsed_args["SEED"])
        print(sample_price_paths[:, 1])
    else
        return rdpp(parsed_args["SEED"])
    end


end

if "" != PROGRAM_FILE && realpath(@__FILE__) == realpath(PROGRAM_FILE)
    main()
end
