# contrastive.jl — Sister hard-negatives as a stream view (not a loss).

"""One labeled read plus causal V index and nearest same-gene sister indices."""
struct ContrastiveExample
    read::LabeledRead
    positive_v::Int
    negative_v::Vector{Int}
end

"""
Wrap a [`ReadStream`](@ref) and attach V-allele index + nearest sisters.

`sister_idx` is `(k, n_v)` from [`nearest_sister_indices`](@ref).
"""
mutable struct ContrastiveStream{S<:ReadStream}
    stream::S
    v_alleles::Vector{Allele}
    name_to_idx::Dict{String,Int}
    sister_idx::Matrix{Int}
end

function ContrastiveStream(stream::ReadStream, db::GermlineDB; sister_k::Integer = 4)
    sis = nearest_sister_indices(db.v; k = Int(sister_k))
    name_to_idx = Dict{String,Int}(a.name => i for (i, a) in enumerate(db.v))
    ContrastiveStream(stream, db.v, name_to_idx, sis)
end

ContrastiveStream(gen::ReadGenerator; sister_k::Integer = 4, seed::Integer = 42) =
    ContrastiveStream(ReadStream(gen; seed), gen.db; sister_k)

Base.IteratorSize(::Type{<:ContrastiveStream}) = Base.IsInfinite()
Base.eltype(::Type{<:ContrastiveStream}) = ContrastiveExample

function Base.iterate(cs::ContrastiveStream, ::Nothing = nothing)
    r = cs.stream.gen(cs.stream.rng)
    pos = get(cs.name_to_idx, r.v_call, 0)
    negs = Int[]
    if pos > 0
        k = size(cs.sister_idx, 1)
        @inbounds for i in 1:k
            j = cs.sister_idx[i, pos]
            j > 0 && push!(negs, j)
        end
    end
    ContrastiveExample(r, pos, negs), nothing
end

function take_batch(cs::ContrastiveStream, n::Integer)
    n = Int(n)
    out = Vector{ContrastiveExample}(undef, n)
    @inbounds for i in 1:n
        out[i] = first(iterate(cs))
    end
    out
end
