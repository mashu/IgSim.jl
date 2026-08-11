# noise.jl — Substitutions + rare indels on VDJ body only (flanks untouched).

"""True for V/D/J/N labels (not flank PAD)."""
is_body_label(lab::Integer) = UInt8(lab) != UInt8(REG_PAD)

"""
Apply body-only substitutions and sparse indels.

Returns `(sequence, labels, n_subs)`. Indels may shift labels; inserts copy
neighbor body label.
"""
function apply_body_noise(rng::AbstractRNG, seq::AbstractString,
                          labels::Vector{UInt8},
                          sub_rate::Real, indel_rate::Real)
    chars = collect(seq)
    labs = copy(labels)
    length(chars) == length(labs) || throw(DimensionMismatch("seq/labels"))
    # indels first (walk copy)
    if Float64(indel_rate) > 0
        out_c = Char[]
        out_l = UInt8[]
        sizehint!(out_c, length(chars) + 8)
        @inbounds for i in eachindex(chars)
            if is_body_label(labs[i]) && rand(rng) < Float64(indel_rate)
                if rand(rng) < 0.5 && length(out_c) > 10
                    # deletion: skip
                    continue
                else
                    push!(out_c, chars[i])
                    push!(out_l, labs[i])
                    push!(out_c, rand(rng, DNA_BASES))
                    push!(out_l, labs[i])
                    continue
                end
            end
            push!(out_c, chars[i])
            push!(out_l, labs[i])
        end
        chars = out_c
        labs = out_l
    end
    n_subs = 0
    r = Float64(sub_rate)
    if r > 0
        @inbounds for i in eachindex(chars)
            if is_body_label(labs[i]) && rand(rng) < r
                b = chars[i]
                choices = b == 'A' ? ('C', 'G', 'T') :
                          b == 'C' ? ('A', 'G', 'T') :
                          b == 'G' ? ('A', 'C', 'T') : ('A', 'C', 'G')
                chars[i] = rand(rng, choices)
                n_subs += 1
            end
        end
    end
    String(chars), labs, n_subs
end
