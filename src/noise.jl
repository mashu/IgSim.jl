# noise.jl — SHM + Illumina substitutions on VDJ body only (flanks untouched).

"""True for V/D/J/N labels (not flank PAD)."""
is_body_label(lab::Integer) = UInt8(lab) != UInt8(REG_PAD)

function substitute_base(rng::AbstractRNG, b::Char)
    b == 'A' ? rand(rng, ('C', 'G', 'T')) :
    b == 'C' ? rand(rng, ('A', 'G', 'T')) :
    b == 'G' ? rand(rng, ('A', 'C', 'T')) : rand(rng, ('A', 'C', 'G'))
end

"""
Apply body-only SHM substitutions and sparse indels.

Returns `(sequence, labels, n_subs)`. Indels may shift labels; inserts copy
neighbor body label. Illumina error is a separate later pass.
"""
function apply_body_noise(rng::AbstractRNG, seq::AbstractString,
                          labels::Vector{UInt8},
                          sub_rate::Real, indel_rate::Real)
    chars = collect(seq)
    labs = copy(labels)
    length(chars) == length(labs) || throw(DimensionMismatch("seq/labels"))
    if Float64(indel_rate) > 0
        out_c = Char[]
        out_l = UInt8[]
        sizehint!(out_c, length(chars) + 8)
        @inbounds for i in eachindex(chars)
            if is_body_label(labs[i]) && rand(rng) < Float64(indel_rate)
                if rand(rng) < 0.5 && length(out_c) > 10
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
                chars[i] = substitute_base(rng, chars[i])
                n_subs += 1
            end
        end
    end
    String(chars), labs, n_subs
end

"""
Apply Illumina-like substitutions on the VDJ body after SHM.

Rate interpolates linearly from `err.p5` at the 5′ body base to `err.p3` at
the 3′ body base. Flanks are left unchanged. Returns `(sequence, n_subs)`.
Labels are not shifted (substitutions only).
"""
function apply_illumina_error(rng::AbstractRNG, seq::AbstractString,
                              labels::Vector{UInt8}, err::IlluminaError)
    p5, p3 = err.p5, err.p3
    (p5 <= 0 && p3 <= 0) && return String(seq), 0
    chars = collect(seq)
    length(chars) == length(labels) || throw(DimensionMismatch("seq/labels"))
    n_body = 0
    @inbounds for lab in labels
        is_body_label(lab) && (n_body += 1)
    end
    n_body == 0 && return String(seq), 0
    n_err = 0
    k = 0
    @inbounds for i in eachindex(chars)
        is_body_label(labels[i]) || continue
        k += 1
        t = n_body == 1 ? 0.0 : (k - 1) / (n_body - 1)
        p = p5 + t * (p3 - p5)
        if p > 0 && rand(rng) < p
            chars[i] = substitute_base(rng, chars[i])
            n_err += 1
        end
    end
    String(chars), n_err
end
