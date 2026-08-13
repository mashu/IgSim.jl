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
    shm(p, lo, hi) = GatedRate(p, Uniform(lo, hi))

Convenience for a zero-inflated Uniform SHM prior. Equivalent to constructing
[`GatedRate`](@ref) yourself; `lo == hi` uses `Dirac(lo)` instead of Uniform.

- `p` — fraction of reads that receive any SHM (the rest have rate `0`)
- `lo`, `hi` — per-base substitution rate on the VDJ body when SHM is on
"""
function shm(p::Real, lo::Real, hi::Real)
    lo64, hi64 = Float64(lo), Float64(hi)
    lo64 > hi64 && throw(ArgumentError("shm lo must be ≤ hi, got $lo64 > $hi64"))
    rate = lo64 == hi64 ? Dirac(lo64) : Uniform(lo64, hi64)
    GatedRate(p, rate)
end

"""
    shm_none() = GatedRate(0.0, Dirac(0.0))

Unmutated / naive IgM: no substitutions.
"""
shm_none() = GatedRate(0.0, Dirac(0.0))

"""
    shm_igm() = GatedRate(0.40, Uniform(0.005, 0.07))

Typical IgM-sorted repertoire: ~40% of reads with biological SHM, per-base rate
`Uniform(0.5%, 7%)` when on (library-wide mean ≈ 1.5%). Single-mismatch noise
is left to the Illumina layer.
"""
shm_igm() = GatedRate(0.40, Uniform(0.005, 0.07))

"""
    shm_igg() = GatedRate(0.97, Uniform(0.015, 0.13))

Switched / memory IgG: ~97% mutated, per-base rate `Uniform(1.5%, 13%)` when on
(library-wide mean ≈ 7%).
"""
shm_igg() = GatedRate(0.97, Uniform(0.015, 0.13))

"""
    IlluminaError(p5, p3 = p5)

Always-on Illumina substitutions on the VDJ body **after** SHM. Rate ramps
linearly from `p5` at the 5′ body nucleotide to `p3` at the 3′ body nucleotide.
Pass a single value for a uniform rate (the merged paired-end case). Flanks
are skipped — they are already random DNA.

Typical Q30 ≈ `0.001`. Illumina indels are omitted (orders of magnitude rarer
than substitutions). Use [`illumina_miseq`](@ref) / [`illumina_none`](@ref).
"""
struct IlluminaError
    p5::Float64
    p3::Float64
    function IlluminaError(p5::Real, p3::Real = p5)
        a, b = Float64(p5), Float64(p3)
        (0 <= a <= 1 && 0 <= b <= 1) ||
            throw(ArgumentError("IlluminaError rates must be in [0, 1], got $a, $b"))
        new(a, b)
    end
end

"""No instrument substitutions: `IlluminaError(0.0)`."""
illumina_none() = IlluminaError(0.0)

"""
    illumina_miseq() = IlluminaError(0.001)

Merged MiSeq paired-end, uniform ~Q30 (`0.1%`) on the VDJ body. Overlap merge
(FLASH/pRESTO-style, keep higher-Q base) cancels the single-end 3′ decay, so
this is not a 5′→3′ ramp. For unmerged single-end use
`IlluminaError(0.001, 0.004)`.
"""
illumina_miseq() = IlluminaError(0.001)

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
- `illumina_error` — always-on instrument substitutions ([`IlluminaError`](@ref))
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
    illumina_error::IlluminaError
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

Package defaults: empirical trim/N/D-inclusion margins, **IgM-like** gated SHM
([`shm_igm`](@ref)), MiSeq-like Illumina error ([`illumina_miseq`](@ref)), and
domain-randomized flanks. Use [`shm_igg`](@ref) for switched memory, or
[`shm`](@ref) / [`GatedRate`](@ref) to set SHM `p` and the rate range yourself.

| Keyword | Default | Role |
|:--------|:--------|:-----|
| `v_trim_5p` | `DiscreteUniform(0, 0)` | V 5′ trim (nt) |
| `v_trim_3p` | `DiscreteUniform(0, 10)` | V 3′ trim (nt) |
| `d_trim_5p` | `Geometric(0.20)` | D 5′ trim (nt; mean = 4, unbounded tail) |
| `d_trim_3p` | `Geometric(0.21)` | D 3′ trim (nt; mean ≈ 3.8, unbounded tail) |
| `j_trim_5p` | `DiscreteUniform(0, 14)` | J 5′ trim (nt) |
| `n1_length` | `DiscreteUniform(0, 18)` | N1 addition length (nt) |
| `n2_length` | `DiscreteUniform(0, 18)` | N2 addition length (nt) |
| `include_d` | `Bernoulli(0.99)` | rare VJ-only skip |
| `flank_5p` | [`default_flank_5p`](@ref) | 5′ unread flank |
| `flank_3p` | [`default_flank_3p`](@ref) | 3′ unread flank |
| `body_error_rate` | `GatedRate(0.40, Uniform(0.005, 0.07))` | SHM ([`shm_igm`](@ref)) |
| `indel_rate` | `GatedRate(0.08, Uniform(0.0, 0.002))` | indel gate + rate |
| `illumina_error` | `IlluminaError(0.001)` | instrument ([`illumina_miseq`](@ref)) |
| `min_length` | `80` | minimum accepted read length (nt) |
| `max_length` | `900` | maximum accepted read length (nt) |
| `max_retries` | `32` | recombine attempts before error |

```julia
params = train_params(; body_error_rate = GatedRate(0.97, Uniform(0.015, 0.13)))
gen = ReadGenerator(db; params)
```
"""
function train_params(;
                      v_trim_5p = DiscreteUniform(0, 0),
                      v_trim_3p = DiscreteUniform(0, 10),
                      d_trim_5p = Geometric(0.20),
                      d_trim_3p = Geometric(0.21),
                      j_trim_5p = DiscreteUniform(0, 14),
                      n1_length = DiscreteUniform(0, 18),
                      n2_length = DiscreteUniform(0, 18),
                      include_d = Bernoulli(0.99),
                      flank_5p = default_flank_5p(),
                      flank_3p = default_flank_3p(),
                      body_error_rate = shm_igm(),
                      indel_rate = GatedRate(0.08, Uniform(0.0, 0.002)),
                      illumina_error = illumina_miseq(),
                      min_length::Integer = 80,
                      max_length::Integer = 900,
                      max_retries::Integer = 32)
    SimParams(
        v_trim_5p, v_trim_3p, d_trim_5p, d_trim_3p, j_trim_5p,
        n1_length, n2_length, include_d, flank_5p, flank_3p,
        body_error_rate, indel_rate, illumina_error,
        Int(min_length), Int(max_length), Int(max_retries),
    )
end

"""Easy: light trims, no SHM (unmutated / naive IgM)."""
easy_params() = train_params(
    v_trim_5p = DiscreteUniform(0, 0),
    v_trim_3p = DiscreteUniform(0, 4),
    d_trim_5p = DiscreteUniform(0, 2),
    d_trim_3p = DiscreteUniform(0, 2),
    j_trim_5p = DiscreteUniform(0, 3),
    n1_length = DiscreteUniform(0, 4),
    n2_length = DiscreteUniform(0, 4),
    include_d = Bernoulli(0.99),
    body_error_rate = shm_none(),
    indel_rate = GatedRate(0.0, Dirac(0.0)),
    flank_5p = DomainFlankMix(DiscreteUniform(10, 55),
                              DiscreteUniform(70, 130), 0.5),
    flank_3p = DomainFlankMix(DiscreteUniform(40, 100),
                              DiscreteUniform(80, 140), 0.5),
)

"""Mid: IgM-library SHM (same as [`train_params`](@ref) defaults)."""
mid_params() = train_params()

"""
Hard: heavier trim + IgG-level gated SHM ([`shm_igg`](@ref)).

~97% of hard reads mutated; rate ~1.5–13% when on.
"""
hard_params() = train_params(
    v_trim_5p = DiscreteUniform(0, 1),
    v_trim_3p = DiscreteUniform(0, 18),
    j_trim_5p = DiscreteUniform(0, 14),
    n1_length = NegativeBinomial(3, 0.35),
    n2_length = NegativeBinomial(3, 0.35),
    body_error_rate = shm_igg(),
    indel_rate = GatedRate(0.15, Uniform(0.0, 0.004)),
    flank_5p = DomainFlankMix(DiscreteUniform(10, 60),
                              DiscreteUniform(70, 150), 0.6),
    flank_3p = DomainFlankMix(DiscreteUniform(35, 110),
                              DiscreteUniform(80, 160), 0.6),
)

"""Assay-matched short flanks (robustness / ablation); IgM-library SHM."""
short_flank_params() = train_params(
    flank_5p = DiscreteUniform(15, 55),
    flank_3p = DiscreteUniform(40, 100),
    max_length = 550,
)
