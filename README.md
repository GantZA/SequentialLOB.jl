# A Sequential Model of the Latent-Order Book Implemented in Julia

## Authors
* Michael Gant
* Tim Gebbie

## Acknowledgements

We would like to thank Byron Jacobs and Chris Angstmann for their assistance and advice with the formulation of a numerical solution of the SPDE necessary for our implementation. We would also like to thank Donovan Platt for his research into calibrating the latent order book model and the code patterns for the numerical solution and calibration techniques. Michael Gant would like to thank Tim Gebbie for his supervision which entailed numerous discussions, valuable advice and resources.


## Arguments

The `SLOB` object is created by instantiating the struct with the relevant parameters. Once instantiated, the object can be called as a function, with a seed as the only argument and mid-price paths will be renturned.

There are 2 command line interfaces which can be accesed. The commands and arguments are detailed below.
1) Running `$ julia src/main.jl arg1 arg2 ...`
2) After building the executable, `slob arg1 arg2 ...`

Each command line interface uses the same positional arguments which are:
* SEED :: Integer - The seed used for any random number generation. This ensures that paths are reproducible. A value of -1 will generate and use random seeds.
* num_paths :: Integer - The number of price paths to simulate
* T :: Integer - Number of time periods that are simulated
* p₀ :: Float - The initial mid-price that each simulation begins at
* M :: Integer - The number of discretized price points used to solve the PDE
* L :: Float - The length of the price grid. L and p₀ determine x₀ and xₘ
* D :: Float - The Diffusion coefficient in the PDE
* σ :: Float - The scaling value for the Stochastic Drift term
* nu :: Float - The latent order cancellation rate, should be set to 0.0 in the Simple LOB model
* α :: Float - The scale parameter (equal to mean) for the Exponential waiting time between Latent Order Book denisty resets
* λ :: Float - Source Term function parameter 1
* μ :: Float - Source Term function parameter 2
## Example Usage

### Julia Terminal

```
julia> using SequentialLOB
julia> slob = SLOB(num_paths=1 ,T=2300, p₀=238.745,
  M=200, L=100.0, D=1.0, σ=0.1, nu=0.1, α=10.0, λ=1.0, μ=0.2)
julia> slob(45)

```

### Shell
```
$ julia src/main.jl 45 1 2300 238.745 200 100.0 1.0 0.1 0.1 10.0 1.0 0.2
```
