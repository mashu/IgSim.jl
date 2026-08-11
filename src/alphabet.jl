# alphabet.jl — Minimal DNA alphabet.

const DNA_BASES = ('A', 'C', 'G', 'T')

struct NucleotideAlphabet end

vocab_size(::NucleotideAlphabet) = 5  # PAD + ACGT
pad_token(::NucleotideAlphabet) = Int32(1)

function tokenize(::NucleotideAlphabet, seq::AbstractString)
    out = Vector{Int32}(undef, length(seq))
    @inbounds for (i, c) in enumerate(seq)
        out[i] = c == 'A' || c == 'a' ? Int32(2) :
                 c == 'C' || c == 'c' ? Int32(3) :
                 c == 'G' || c == 'g' ? Int32(4) :
                 c == 'T' || c == 't' ? Int32(5) : Int32(1)
    end
    out
end

function detokenize(::NucleotideAlphabet, tokens::AbstractVector{<:Integer})
    chars = Vector{Char}(undef, length(tokens))
    @inbounds for i in eachindex(tokens)
        t = Int(tokens[i])
        chars[i] = t == 2 ? 'A' : t == 3 ? 'C' : t == 4 ? 'G' : t == 5 ? 'T' : 'N'
    end
    String(chars)
end

random_dna(rng::AbstractRNG, n::Integer) =
    String([rand(rng, DNA_BASES) for _ in 1:Int(n)])
