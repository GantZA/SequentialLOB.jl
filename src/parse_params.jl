using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "SEED"
            help = "Seed for randomness"
            arg_type = Int
            default = 1
        "num_paths"
            help = "Number of price paths to generate"
            arg_type = Int
            default = 1
        "T"
            help = "Number of time periods"
            arg_type = Int
            default = 100
        "p₀"
            arg_type = Float64
            default = 100.0
        "M"
            help = "Number of price points"
            arg_type = Int
            default = 100
        "β"
            arg_type = Float64
            default = 1.0
        "L"
            arg_type = Float64
            default = 50.0
        "D"
            arg_type = Float64
            default = 4.0
        "σ"
            arg_type = Float64
            default = 1.0
        "nu"
            arg_type = Float64
            default = 0.0
        "α"
            arg_type = Float64
            default = 1.0
        "λ"
            arg_type = Float64
            default = 1.0
        "μ"
            arg_type = Float64
            default = 0.5
    end

    return parse_args(s)
end
