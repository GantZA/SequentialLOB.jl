using SequentialLOB
using Test

include("../src/main.jl")

@testset "All Tests" begin
    slob_model = SLOB()
    @testset "Reproducibility" begin
        slob_model_1 = SLOB(1, 200, 238.745, 500, 100.0, 2.0,
            1.0, 0.1, 20.0, SourceTerm(1.0, 0.5))
        slob_model_2 = SLOB(1, 200, 238.745, 500, 100.0, 2.0,
            1.0, 0.1, 20.0, SourceTerm(1.0, 0.5))
        @test all(slob_model_1(45) .== slob_model_2(45))

        @test all(slob_model_1(51) .== slob_model_2(51))

        @test all(slob_model_1(4565756745) .== slob_model_2(4565756745))
    end;

    @testset "Default Values" begin
        lob_densities, price_paths, mid_price_bars, P⁺s, P⁻s = slob_model()
        @test size(mid_price_bars) == (101,1)
    end;

    @testset "Command Line Parse Default Arguments" begin
        @test parse_commandline() == Dict("μ" => 0.5,
            "num_paths" => 1,
            "T" => 100,
            "p₀" => 100.0,
            "λ" => 1.0,
            "SEED" => 1,
            "nu" => 0.1,
            "σ" => 1.0,
            "D" => 4.0,
            "M" => 100,
            "L" => 50.0,
            "α" => 20.0)
    end

    @testset "Main Default Arguments" begin
        @test all(main("return") .== slob_model(1))
    end
end
