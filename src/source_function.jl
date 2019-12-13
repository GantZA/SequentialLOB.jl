# """
# Source Term Function
# ```julia
# source_term = SourceTerm(λ, μ)
# source_term(x, p)
# ```
# where x is a discritzed price point, source_params is a vector of parameters for
# the source term and p is the last "observed" mid price point
#
# The function returns a single scalar output
#
# """

mutable struct SourceTerm
    λ::Float64
    μ::Float64
end


# function (st::SourceTerm)(x::Float64,p::Float64)
#     return st.λ*(p-x)*exp(-st.μ*(p-x)^2)
# end

function (st::SourceTerm)(x::Float64,p::Float64)
    return st.λ*tanh(st.μ*(p-x))
end
