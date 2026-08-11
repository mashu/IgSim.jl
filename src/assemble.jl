# assemble.jl — Concatenate parts + flanks into sequence and labels.

function append_seg!(seq::Vector{Char}, labels::Vector{UInt8},
                     piece::AbstractString, lab::RegionLabel)
    isempty(piece) && return
    for c in piece
        push!(seq, c)
        push!(labels, UInt8(lab))
    end
end

"""
Build full read: flank5 + V + N1 + [D] + N2 + J + flank3.

Returns `(sequence, labels, flank5_len, flank3_len)`.
"""
function assemble(rng::AbstractRNG, parts::RecomboParts, params::SimParams)
    f5 = Int(rand(rng, params.flank_5p))
    f3 = Int(rand(rng, params.flank_3p))
    seq = Char[]
    labels = UInt8[]
    sizehint!(seq, 600)
    sizehint!(labels, 600)
    append_seg!(seq, labels, random_dna(rng, f5), REG_PAD)
    append_seg!(seq, labels, parts.v_body, REG_V)
    append_seg!(seq, labels, parts.n1, REG_N1)
    if parts.has_d
        append_seg!(seq, labels, parts.d_body, REG_D)
        append_seg!(seq, labels, parts.n2, REG_N2)
    end
    append_seg!(seq, labels, parts.j_body, REG_J)
    append_seg!(seq, labels, random_dna(rng, f3), REG_PAD)
    String(seq), labels, f5, f3
end
