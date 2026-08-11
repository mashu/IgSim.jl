# IgSim.jl

[![Build Status](https://github.com/mashu/IgSim.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mashu/IgSim.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mashu/IgSim.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/IgSim.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mashu.github.io/IgSim.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/IgSim.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Streaming V(D)J read simulator for training: uniform allele sampling, labeled
spans, holdout splits, and contrastive sister negatives.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/mashu/IgSim.jl")
```

## Usage

```julia
using IgSim, Random

db = load_germline(;
    v = "path/to/V.fasta",
    d = "path/to/D.fasta",   # optional: omit or `nothing` for VJ-only
    j = "path/to/J.fasta",
)

gen = ReadGenerator(db)                      # train_params() defaults
reads = take_batch(gen, 64; seed = 1)
batch = collate_reads(reads)

# contrastive view: causal V index + nearest same-gene sisters
cs = ContrastiveStream(gen; sister_k = 4)
examples = take_batch(cs, 64)
```

### Holdout (open-set)

```julia
train_db, held = holdout_alleles(full_db, held_v_names, String[], String[])
train_gen = ReadGenerator(train_db)                 # never samples held V
held_gen  = HoldoutVGenerator(full_db, held_v_names) # gold V ∈ held set
# evaluation gallery = full_db (or merge_germline(train_db, held))
```

### Custom parameters

Defaults come from `train_params`. Override any keyword (full list and defaults
in the [Parameters](https://mashu.github.io/IgSim.jl/dev/parameters/) docs):

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
| `flank_5p` | `DomainFlankMix(...)` | 5′ unread flank |
| `flank_3p` | `DomainFlankMix(...)` | 3′ unread flank |
| `body_error_rate` | `Uniform(0.0, 0.06)` | substitution rate on VDJ body |
| `indel_rate` | `Uniform(0.0, 0.002)` | indel rate on VDJ body |
| `min_length` | `80` | min accepted length (nt) |
| `max_length` | `900` | max accepted length (nt) |
| `max_retries` | `32` | recombine attempts before error |

```julia
using Distributions
params = train_params(;
    flank_5p = DomainFlankMix(DiscreteUniform(5, 40), DiscreteUniform(40, 120), 0.6),
    body_error_rate = Uniform(0.0, 0.03),
)
gen = ReadGenerator(db; params)
```

## License

MIT — see [LICENSE](LICENSE).

