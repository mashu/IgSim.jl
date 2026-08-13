# recombine.jl — DJ first, then V+DJ (uniform allele draws).

"""Parts of a recombined read before flank/noise (ASCII DNA strings)."""
struct RecomboParts
    v_body::String
    n1::String
    d_body::String
    n2::String
    j_body::String
    v_call::String
    d_call::Union{String,Missing}
    j_call::String
    has_d::Bool
    v_trim_5::Int
    v_trim_3::Int
    d_trim_5::Int
    d_trim_3::Int
    j_trim_5::Int
end

"""
Sample alleles uniformly and assemble DJ then V+DJ with trims/N.

Order: sample J → (optional D) → trim D3′/J5′ → N2 → DJ;
then sample V → trim V3′ / DJ5′(via D5 or leftover) → N1 → VDJ.
"""
function recombine(rng::AbstractRNG, db::GermlineDB, params::SimParams)
    j_all = sample_allele(rng, db.j)
    j_trim = Int(rand(rng, params.j_trim_5p))
    j_body = trim_5p_ascii(j_all.sequence, j_trim)

    has_d = !isempty(db.d) && rand(rng, params.include_d)
    d_call = missing
    d_body = ""
    d_trim_5 = 0
    d_trim_3 = 0
    n2 = ""
    if has_d
        d_all = sample_allele(rng, db.d)
        d_trim_5 = Int(rand(rng, params.d_trim_5p))
        d_trim_3 = Int(rand(rng, params.d_trim_3p))
        d_body = trim_5p_ascii(trim_3p_ascii(d_all.sequence, d_trim_3), d_trim_5)
        # Keep the sampled D label whenever any nucleotide remains, including
        # 1–2 nt remnants. Call no-D only when trimming ate the whole D.
        if isempty(d_body)
            has_d = false
        else
            d_call = d_all.name
            n2 = sample_n(rng, Int(rand(rng, params.n2_length)))
        end
    end

    v_all = sample_allele(rng, db.v)
    v_trim_5 = Int(rand(rng, params.v_trim_5p))
    v_trim_3 = Int(rand(rng, params.v_trim_3p))
    v_body = trim_5p_ascii(trim_3p_ascii(v_all.sequence, v_trim_3), v_trim_5)
    isempty(v_body) && return nothing
    isempty(j_body) && return nothing

    n1 = sample_n(rng, Int(rand(rng, params.n1_length)))

    RecomboParts(v_body, n1, d_body, n2, j_body,
                 v_all.name, d_call, j_all.name, has_d,
                 v_trim_5, v_trim_3, d_trim_5, d_trim_3, j_trim)
end
