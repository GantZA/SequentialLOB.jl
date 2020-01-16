mutable struct SourceTerm
    λ::Float64
    μ::Float64
end


function (st::SourceTerm)(x::Float64,p::Float64)
    return st.λ*tanh(st.μ*(p-x))
end
