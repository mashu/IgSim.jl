# Parameters

All simulation knobs live in [`SimParams`](@ref), built by
[`train_params`](@ref). Pass keyword overrides to `train_params`, then
`ReadGenerator(db; params)`.

```@docs
train_params
SimParams
DomainFlankMix
GatedRate
shm
shm_none
shm_igm
shm_igg
IlluminaError
illumina_none
illumina_miseq
default_flank_5p
default_flank_3p
```

## Override checklist

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

Flanks use [`DomainFlankMix`](@ref): with probability `p_short` sample from
`short`, otherwise from `long`.

A read is labeled no-D only if `include_d` skips the segment or exonuclease
leaves **zero** D nucleotides. Any remnant of 1 nt or more keeps the originally
sampled D label. D 5′/3′ deletions are `Geometric` (mean ≈ 4 nt, unbounded
tail) — the **same** process on every D gene, then clamped by available
length. Short Ds are fully eaten more often as a consequence, not via a
special case. This is causal gold, not an IgBLAST minimum consecutive-match
cutoff. The mean is set so IgBLAST / SwiftIG empty-`d_call` rates on sim
land near the same tools on real reads; gold emptiness stays higher than
either caller reports.

## Error model

SHM and Illumina both apply to the **VDJ body** (not flanks). SHM first, then
Illumina. Presets are zero-arg aliases — the implementation is the
`GatedRate` / `IlluminaError` on the right.

**SHM** is [`GatedRate`](@ref): with probability `p` draw a per-base rate from
`rate`, otherwise `0`. Each body base then substitutes independently.

| Alias | Implementation | Typical use |
|:------|:---------------|:------------|
| [`shm_none`](@ref) | `GatedRate(0.0, Dirac(0.0))` | naive / unmutated IgM |
| [`shm_igm`](@ref) | `GatedRate(0.40, Uniform(0.005, 0.07))` | IgM-sorted library (default) |
| [`shm_igg`](@ref) | `GatedRate(0.97, Uniform(0.015, 0.13))` | switched memory IgG |
| [`shm`](@ref) `(p, lo, hi)` | `GatedRate(p, Uniform(lo, hi))` | custom mix |

This is a training prior (naive peak + mutated clones), not a fit to one
repertoire. Uniform over mutation depth is a flat prior; for a heavy tail use
`GatedRate(p, LogNormal(μ, σ))`. No AID hotspot / CDR targeting.

**Illumina** ([`IlluminaError`](@ref)) is always-on and applied after SHM.
Default `IlluminaError(0.001)` is uniform ~Q30 on a **merged paired-end**
contig (overlap keeps the higher-Q base, so single-mate errors are usually
dropped). For unmerged single-end use a ramp, e.g. `IlluminaError(0.001, 0.004)`.
Disable with `IlluminaError(0.0)` / [`illumina_none`](@ref). Substitutions only.

Replace either flank with a single `DiscreteUniform` (or any sampleable) if you
want a simpler length model.
