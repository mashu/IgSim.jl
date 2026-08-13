# recombine.jl — DJ first, optional D+DJ (tandem DD), then V+(D)DJ.

"""Parts of a recombined read before flank/noise (ASCII DNA strings)."""
struct RecomboParts
    v_body::String
    n1::String
    d_body::String
    n3::String
    d2_body::String
    n2::String
    j_body::String
    v_call::String
    d_call::Union{String,Missing}
    j_call::String
    has_d::Bool
    multi_d::Bool
    v_trim_5::Int
    v_trim_3::Int
    d_trim_5::Int
    d_trim_3::Int
    d2_trim_5::Int
    d2_trim_3::Int
    j_trim_5::Int
end

"""Independent 5′/3′ exonuclease on one D allele."""
function trim_d(rng::AbstractRNG, params::SimParams, seq::AbstractString)
    t5 = Int(rand(rng, params.d_trim_5p))
    t3 = Int(rand(rng, params.d_trim_3p))
    body = trim_5p_ascii(trim_3p_ascii(seq, t3), t5)
    body, t5, t3
end

"""
Sample alleles uniformly and assemble DJ, optional DDJ, then V+(D)DJ.

Order: sample J → D → trim D3′/J5′ → N2 → DJ;
with small probability a second D joins that DJ (trim both ends, N3) → DDJ;
then sample V → trim V3′ / (D)DJ5′ → N1 → V(D)DJ.

Gold `d_call` is a single allele, or a comma-separated 5′ then 3′ pair
(same string form as assignment-tool multi-calls). `multi_d` is true only
when both D remnants are non-empty.
"""
function recombine(rng::AbstractRNG, db::GermlineDB, params::SimParams)
    j_all = sample_allele(rng, db.j)
    j_trim = Int(rand(rng, params.j_trim_5p))
    j_body = trim_5p_ascii(j_all.sequence, j_trim)

    has_d = !isempty(db.d) && rand(rng, params.include_d)
    d_call = missing
    d_body = ""
    d2_body = ""
    d_trim_5 = 0
    d_trim_3 = 0
    d2_trim_5 = 0
    d2_trim_3 = 0
    n2 = ""
    n3 = ""
    multi_d = false
    if has_d
        d3_i = rand(rng, 1:length(db.d))
        d3_all = db.d[d3_i]
        d3_body, d3_t5, d3_t3 = trim_d(rng, params, d3_all.sequence)
        # Keep the sampled D label whenever any nucleotide remains, including
        # 1–2 nt remnants. Call no-D only when trimming ate the whole D.
        if isempty(d3_body)
            has_d = false
        else
            n2 = sample_n(rng, Int(rand(rng, params.n2_length)))
            d_body = d3_body
            d_trim_5, d_trim_3 = d3_t5, d3_t3
            d_call = d3_all.name
            want_dd = length(db.d) >= 2 && rand(rng, params.include_dd)
            if want_dd
                _, d5_all = sample_other_allele(rng, db.d, d3_i)
                d5_body, d5_t5, d5_t3 = trim_d(rng, params, d5_all.sequence)
                if !isempty(d5_body)
                    # Second recombination: D + DJ → DDJ. 5′ D is the new one.
                    n3 = sample_n(rng, Int(rand(rng, params.n2_length)))
                    d2_body = d3_body
                    d2_trim_5, d2_trim_3 = d3_t5, d3_t3
                    d_body = d5_body
                    d_trim_5, d_trim_3 = d5_t5, d5_t3
                    d_call = d5_all.name * "," * d3_all.name
                    multi_d = true
                end
            end
        end
    end

    v_all = sample_allele(rng, db.v)
    v_trim_5 = Int(rand(rng, params.v_trim_5p))
    v_trim_3 = Int(rand(rng, params.v_trim_3p))
    v_body = trim_5p_ascii(trim_3p_ascii(v_all.sequence, v_trim_3), v_trim_5)
    isempty(v_body) && return nothing
    isempty(j_body) && return nothing

    n1 = sample_n(rng, Int(rand(rng, params.n1_length)))

    RecomboParts(v_body, n1, d_body, n3, d2_body, n2, j_body,
                 v_all.name, d_call, j_all.name, has_d, multi_d,
                 v_trim_5, v_trim_3, d_trim_5, d_trim_3,
                 d2_trim_5, d2_trim_3, j_trim)
end
