# read.jl — LabeledRead + ReadGenerator.

"""One simulated read with causal allele labels and region path.

`d_call` is a single allele or a comma-separated 5′ then 3′ pair (no spaces),
the same string form assignment tools use for multiple D calls. `multi_d` is
`true` only for tandem DD (both remnants non-empty), so a comma in `d_call`
is gold double-D rather than annotator ambiguity.
"""
struct LabeledRead
    sequence::String
    labels::Vector{UInt8}
    spans::NamedTuple
    v_call::String
    d_call::Union{String,Missing}
    j_call::String
    has_d::Bool
    multi_d::Bool
    flank5::Int
    flank3::Int
    v_trim_5::Int
    v_trim_3::Int
    d_trim_5::Int
    d_trim_3::Int
    d2_trim_5::Int
    d2_trim_3::Int
    j_trim_5::Int
    n_errors::Int
    n_illumina::Int
    body_error_rate::Float64
end

"""Callable generator: `gen(rng) -> LabeledRead`."""
struct ReadGenerator{P<:SimParams}
    db::GermlineDB
    params::P
end

ReadGenerator(db::GermlineDB; params::SimParams = train_params()) =
    ReadGenerator(db, params)

function (gen::ReadGenerator)(rng::AbstractRNG)
    params = gen.params
    for _ in 1:params.max_retries
        parts = recombine(rng, gen.db, params)
        isnothing(parts) && continue
        seq0, labels0, f5, f3 = assemble(rng, parts, params)
        L = length(seq0)
        (L < params.min_length || L > params.max_length) && continue
        sub_r = Float64(rand(rng, params.body_error_rate))
        ind_r = Float64(rand(rng, params.indel_rate))
        seq, labels, n_err = apply_body_noise(rng, seq0, labels0, sub_r, ind_r)
        seq, n_ill = apply_illumina_error(rng, seq, labels, params.illumina_error)
        spans = spans_from_labels(labels)
        return LabeledRead(seq, labels, spans,
                           parts.v_call, parts.d_call, parts.j_call,
                           parts.has_d, parts.multi_d,
                           f5, f3,
                           parts.v_trim_5, parts.v_trim_3,
                           parts.d_trim_5, parts.d_trim_3,
                           parts.d2_trim_5, parts.d2_trim_3, parts.j_trim_5,
                           n_err, n_ill, sub_r)
    end
    error("ReadGenerator failed after $(params.max_retries) retries")
end

"""
Open-set sim: sample V only from held alleles; D/J from full DB.

`held_v_names` are exact allele strings. Typical protocol: train on
`holdout_alleles(...).train`, evaluate gallery = full DB, sim reads from
`HoldoutVGenerator(full_db, held_names)`.
"""
struct HoldoutVGenerator{G<:ReadGenerator}
    gen::G
    held_v_names::Set{String}
end

function HoldoutVGenerator(db_full::GermlineDB, held_v_names;
                           params::SimParams = train_params())
    held_v = filter_alleles(db_full.v, held_v_names)
    isempty(held_v) && throw(ArgumentError("no held V alleles matched"))
    gen = ReadGenerator(GermlineDB(held_v, db_full.d, db_full.j); params)
    HoldoutVGenerator(gen, Set(String.(a.name for a in held_v)))
end

(h::HoldoutVGenerator)(rng::AbstractRNG) = h.gen(rng)

"""
    MixedReadGenerator(gens, weights)

Per-read mixture over generators (e.g. easy / mid / hard curriculum).
`weights` are relative (need not sum to 1).
"""
struct MixedReadGenerator{G}
    gens::Vector{G}
    weights::Vector{Float64}
    function MixedReadGenerator(gens::Vector{G},
                                weights::AbstractVector{<:Real}) where {G}
        length(gens) == length(weights) || throw(ArgumentError("gens/weights length mismatch"))
        length(gens) >= 1 || throw(ArgumentError("need at least one generator"))
        w = Float64.(weights)
        any(<(0), w) && throw(ArgumentError("weights must be non-negative"))
        sum(w) > 0 || throw(ArgumentError("weights must sum > 0"))
        new{G}(gens, w)
    end
end

"""
Default train mix: naive IgM / IgM-library SHM / IgG tail.

Weights `(0.55, 0.35, 0.10)` → expected fraction of reads with any SHM
≈ `0.55·0 + 0.35·0.40 + 0.10·0.97 ≈ 24%`.
"""
function mixed_curriculum(db::GermlineDB; weights = (0.55, 0.35, 0.10))
    MixedReadGenerator(
        ReadGenerator[
            ReadGenerator(db; params = easy_params()),
            ReadGenerator(db; params = mid_params()),
            ReadGenerator(db; params = hard_params()),
        ],
        collect(Float64, weights),
    )
end

function (gen::MixedReadGenerator)(rng::AbstractRNG)
    total = sum(gen.weights)
    u = rand(rng) * total
    c = 0.0
    @inbounds for i in eachindex(gen.weights)
        c += gen.weights[i]
        u <= c && return gen.gens[i](rng)
    end
    gen.gens[end](rng)
end
