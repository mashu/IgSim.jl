# params.jl — SimParams, gated SHM, and curriculum presets.

"""
Two-component flank-length mixture for unread 5′/3′ contexts.

With probability `p_short`, sample length from `short`; otherwise from `long`.
Used as the default for `flank_5p` / `flank_3p` in [`train_params`](@ref).
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
    GatedRate(p, rate)

Two-level body-noise knob for `body_error_rate` / `indel_rate`:

1. With probability `p`, noise occurs at all.
2. Conditional on that, the per-base rate is drawn from `rate`.

Otherwise the sampled rate is `0` (clean / unmutated read).
"""
struct GatedRate{R}
    p::Float64
    rate::R
    function GatedRate(p::Real, rate)
        (0 <= p <= 1) || throw(ArgumentError("GatedRate p must be in [0, 1], got $p"))
        new{typeof(rate)}(Float64(p), rate)
    end
end

Random.rand(rng::AbstractRNG, g::GatedRate) =
    rand(rng) < g.p ? Float64(rand(rng, g.rate)) : 0.0

"""
    SimParams

Stochastic knobs for recombination, flanks, and body noise. Parametric over
distribution types. Construct with [`train_params`](@ref) or curriculum
helpers [`easy_params`](@ref) / [`mid_params`](@ref) / [`hard_params`](@ref).

# Fields
- `v_trim_5p`, `v_trim_3p` — V end trimming (nt)
- `d_trim_5p`, `d_trim_3p` — D end trimming (nt)
- `j_trim_5p` — J 5′ trimming (nt)
- `n1_length`, `n2_length` — N-addition lengths (nt)
- `include_d` — Bernoulli (or similar) whether the D segment is included
- `flank_5p`, `flank_3p` — unread flank lengths (often [`DomainFlankMix`](@ref))
- `body_error_rate`, `indel_rate` — per-read SHM / indel rates (often [`GatedRate`](@ref))
- `min_length`, `max_length` — accepted assembled length (nt)
- `max_retries` — recombination attempts before error
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
    default_flank_5p() -> DomainFlankMix

Default 5′ unread flank lengths for [`train_params`](@ref): short
`DiscreteUniform(10, 60)`, long `DiscreteUniform(70, 140)`, `p_short = 0.55`.
"""
default_flank_5p() = DomainFlankMix(DiscreteUniform(10, 60),
                                    DiscreteUniform(70, 140), 0.55)

"""
    default_flank_3p() -> DomainFlankMix

Default 3′ unread flank lengths for [`train_params`](@ref): short
`DiscreteUniform(40, 110)`, long `DiscreteUniform(80, 160)`, `p_short = 0.55`.
"""
default_flank_3p() = DomainFlankMix(DiscreteUniform(40, 110),
                                    DiscreteUniform(80, 160), 0.55)

"""
    train_params(; kwargs...) -> SimParams

Package defaults: empirical trim/N/D-inclusion margins, **gated** IgG-like body
noise ([`GatedRate`](@ref)), and domain-randomized flanks.

| Keyword | Default | Role |
|:--------|:--------|:-----|
| `v_trim_5p` | `DiscreteUniform(0, 0)` | V 5′ trim (nt) |
| `v_trim_3p` | `DiscreteUniform(0, 10)` | V 3′ trim (nt) |
| `d_trim_5p` | `DiscreteUniform(0, 12)` | D 5′ trim (nt) |
| `d_trim_3p` | `DiscreteUniform(0, 12)` | D 3′ trim (nt) |
| `j_trim_5p` | `DiscreteUniform(0, 14)` | J 5′ trim (nt) |
| `n1_length` | `DiscreteUniform(0, 18)` | N1 addition length (nt) |
| `n2_length` | `DiscreteUniform(0, 18)` | N2 addition length (nt) |
| `include_d` | `Bernoulli(0.97)` | include D segment |
| `flank_5p` | [`default_flank_5p`](@ref) | 5′ unread flank |
| `flank_3p` | [`default_flank_3p`](@ref) | 3′ unread flank |
| `body_error_rate` | `GatedRate(0.45, Uniform(0.005, 0.04))` | SHM gate + rate |
| `indel_rate` | `GatedRate(0.08, Uniform(0.0, 0.002))` | indel gate + rate |
| `min_length` | `80` | minimum accepted read length (nt) |
| `max_length` | `900` | maximum accepted read length (nt) |
| `max_retries` | `32` | recombine attempts before error |

```julia
params = train_params(; body_error_rate = GatedRate(0.3, Uniform(0.01, 0.05)))
gen = ReadGenerator(db; params)
```
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
                      flank_5p = default_flank_5p(),
                      flank_3p = default_flank_3p(),
                      body_error_rate = GatedRate(0.45, Uniform(0.005, 0.04)),
                      indel_rate = GatedRate(0.08, Uniform(0.0, 0.002)),
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

"""Easy: light trims, no SHM (unmutated / naive-like)."""
easy_params() = train_params(
    v_trim_5p = DiscreteUniform(0, 0),
    v_trim_3p = DiscreteUniform(0, 4),
    d_trim_5p = DiscreteUniform(0, 2),
    d_trim_3p = DiscreteUniform(0, 2),
    j_trim_5p = DiscreteUniform(0, 3),
    n1_length = DiscreteUniform(0, 4),
    n2_length = DiscreteUniform(0, 4),
    include_d = Bernoulli(0.95),
    body_error_rate = GatedRate(0.0, Dirac(0.0)),
    indel_rate = GatedRate(0.0, Dirac(0.0)),
    flank_5p = DomainFlankMix(DiscreteUniform(10, 55),
                              DiscreteUniform(70, 130), 0.5),
    flank_3p = DomainFlankMix(DiscreteUniform(40, 100),
                              DiscreteUniform(80, 140), 0.5),
)

"""Mid: IgG-like mild gated SHM (same as [`train_params`](@ref) defaults)."""
mid_params() = train_params()

"""
Hard: heavier trim + more frequent / deeper gated SHM (affinity-matured tail).

~65% of hard reads mutated; rate ~2–8% when on.
"""
hard_params() = train_params(
    v_trim_5p = DiscreteUniform(0, 1),
    v_trim_3p = DiscreteUniform(0, 18),
    d_trim_5p = DiscreteUniform(0, 12),
    d_trim_3p = DiscreteUniform(0, 12),
    j_trim_5p = DiscreteUniform(0, 14),
    n1_length = NegativeBinomial(3, 0.35),
    n2_length = NegativeBinomial(3, 0.35),
    include_d = Bernoulli(0.85),
    body_error_rate = GatedRate(0.65, Uniform(0.02, 0.08)),
    indel_rate = GatedRate(0.15, Uniform(0.0, 0.004)),
    flank_5p = DomainFlankMix(DiscreteUniform(10, 60),
                              DiscreteUniform(70, 150), 0.6),
    flank_3p = DomainFlankMix(DiscreteUniform(35, 110),
                              DiscreteUniform(80, 160), 0.6),
)

"""Assay-matched short flanks (robustness / ablation); mid gated SHM."""
short_flank_params() = train_params(
    flank_5p = DiscreteUniform(15, 55),
    flank_3p = DiscreteUniform(40, 100),
    max_length = 550,
)
