# params.jl — SimParams and package defaults.

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
    SimParams

Stochastic knobs for recombination, flanks, and body noise. Parametric over
distribution types. Construct with [`train_params`](@ref); every keyword of
`train_params` is a field here.

# Fields
- `v_trim_5p`, `v_trim_3p` — V end trimming (nt)
- `d_trim_5p`, `d_trim_3p` — D end trimming (nt)
- `j_trim_5p` — J 5′ trimming (nt)
- `n1_length`, `n2_length` — N-addition lengths (nt)
- `include_d` — Bernoulli (or similar) whether the D segment is included
- `flank_5p`, `flank_3p` — unread flank lengths (often [`DomainFlankMix`](@ref))
- `body_error_rate`, `indel_rate` — per-read SHM / indel rates on VDJ body
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
    train_params(; kwargs...) -> SimParams

Package defaults: empirical trim/N/D-inclusion margins, modest body-noise
curriculum, and **domain-randomized** flanks (wide short/long mix — not locked
to a single primer/assay unread length).

Pass any of the keywords below to override the default. Values that are sampled
per read should be `Distributions.Sampleable` (or [`DomainFlankMix`](@ref) for
flanks); length bounds are integers.

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
| `flank_5p` | `DomainFlankMix(DiscreteUniform(10, 60), DiscreteUniform(60, 140), 0.5)` | 5′ unread flank |
| `flank_3p` | `DomainFlankMix(DiscreteUniform(20, 80), DiscreteUniform(80, 160), 0.5)` | 3′ unread flank |
| `body_error_rate` | `Uniform(0.0, 0.06)` | substitution rate on VDJ body |
| `indel_rate` | `Uniform(0.0, 0.002)` | indel rate on VDJ body |
| `min_length` | `80` | minimum accepted read length (nt) |
| `max_length` | `900` | maximum accepted read length (nt) |
| `max_retries` | `32` | recombine attempts before error |

```julia
using Distributions
params = train_params(;
    v_trim_3p = DiscreteUniform(0, 5),
    include_d = Bernoulli(1.0),
    body_error_rate = Uniform(0.0, 0.02),
)
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
