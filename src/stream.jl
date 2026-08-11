# stream.jl — Lazy / channel streaming of LabeledRead.

abstract type AbstractReadStream end

"""Infinite iterator of [`LabeledRead`](@ref) from a generator."""
mutable struct ReadStream{G<:ReadGenerator} <: AbstractReadStream
    gen::G
    rng::AbstractRNG
end

ReadStream(gen::ReadGenerator; seed::Integer = 42) =
    ReadStream(gen, MersenneTwister(Int(seed)))

Base.IteratorSize(::Type{<:ReadStream}) = Base.IsInfinite()
Base.eltype(::Type{<:ReadStream}) = LabeledRead

function Base.iterate(s::ReadStream, ::Nothing = nothing)
    s.gen(s.rng), nothing
end

"""Take `n` reads from a stream (or generator+rng)."""
function take_batch(stream::ReadStream, n::Integer)
    n = Int(n)
    out = Vector{LabeledRead}(undef, n)
    @inbounds for i in 1:n
        out[i] = stream.gen(stream.rng)
    end
    out
end

function take_batch(gen::ReadGenerator, n::Integer; seed::Integer = 42)
    take_batch(ReadStream(gen; seed), n)
end

"""
Buffered channel producer (fill in a `@async` / `@spawn` task).

Returns a `Channel{LabeledRead}`. Close the channel to stop.
"""
function channel_stream(gen::ReadGenerator; buffersize::Integer = 256,
                        seed::Integer = 42)
    rng = MersenneTwister(Int(seed))
    Channel{LabeledRead}(Int(buffersize)) do ch
        while true
            put!(ch, gen(rng))
        end
    end
end
