# nadd.jl — Non-templated N segments.

"""Sample a random DNA N-addition of length `L` (empty if L≤0)."""
function sample_n(rng::AbstractRNG, L::Integer)
    L = Int(L)
    L <= 0 && return ""
    random_dna(rng, L)
end
