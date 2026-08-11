# Parameters

All simulation knobs live in [`SimParams`](@ref), built by
[`train_params`](@ref). Pass keyword overrides to `train_params`, then
`ReadGenerator(db; params)`.

```@docs
train_params
SimParams
DomainFlankMix
```

## Override checklist

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

Flanks use [`DomainFlankMix`](@ref): with probability `p_short` sample from
`short`, otherwise from `long`. Replace either flank with a single
`DiscreteUniform` (or any sampleable) if you want a simpler length model.
