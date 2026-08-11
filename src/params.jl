# params.jl — SimParams and package defaults.

"""
Two-component flank-length mixture (short vs longer unread 5′/3′ contexts).
"""
struct DomainFlankMix
    short::DiscreteUniform
    long::DiscreteUniform
    p_short::Float64
end

function Random.rand(rng::AbstractRNG, m::DomainFlankMix)
    rand(rng) < m.p_short ? rand(rng, m.short) : rand(rng, m.long)
end

"""
    SimParams

Stochastic knobs for recombination, flanks, and body noise. Parametric over
distribution types. Construct with [`train_params`](@ref) or override fields.
"""
struct SimParams{VT5,VT3,DT5,DT3,JT5,N1,N2,ID,F5,F3,ER,IR}
    v_trim_5p::VT5
    v_trim_3p::VT3
    d_trim_5p::DT5
    d_trim_3p::DT3
    j_trim_5p::JT5
    n1_length::N1
    n2_length::N2
    include_d::ID
    flank_5p::F5
    flank_3p::F3
    body_error_rate::ER
    indel_rate::IR
    min_length::Int
    max_length::Int
    max_retries::Int
end

"""
    train_params(; kwargs...)

Package defaults: empirical trim/N/D-inclusion margins, modest body-noise
curriculum, and **domain-randomized** flanks (wide short/long mix — not locked
to a single primer/assay unread length). Override any keyword to fine-tune.
"""
function train_params(;
                      v_trim_5p = DiscreteUniform(0, 0),
                      v_trim_3p = DiscreteUniform(0, 10),
                      d_trim_5p = DiscreteUniform(0, 12),
                      d_trim_3p = DiscreteUniform(0, 12),
                      j_trim_5p = DiscreteUniform(0, 14),
                      n1_length = DiscreteUniform(0, 18),
                      n2_length = DiscreteUniform(0, 18),
                      include_d = Bernoulli(0.97),
                      flank_5p = DomainFlankMix(DiscreteUniform(10, 60),
                                                DiscreteUniform(60, 140), 0.5),
                      flank_3p = DomainFlankMix(DiscreteUniform(20, 80),
                                                DiscreteUniform(80, 160), 0.5),
                      body_error_rate = Uniform(0.0, 0.06),
                      indel_rate = Uniform(0.0, 0.002),
                      min_length::Integer = 80,
                      max_length::Integer = 900,
                      max_retries::Integer = 32)
    SimParams(
        v_trim_5p, v_trim_3p, d_trim_5p, d_trim_3p, j_trim_5p,
        n1_length, n2_length, include_d, flank_5p, flank_3p,
        body_error_rate, indel_rate,
        Int(min_length), Int(max_length), Int(max_retries),
    )
end
