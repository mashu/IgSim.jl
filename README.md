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

Defaults come from `train_params`. Override any keyword (full list in the
[Parameters](https://mashu.github.io/IgSim.jl/dev/parameters/) docs). Noise
presets (`shm_igm()`, `illumina_miseq()`, …) are zero-arg aliases for the
distributions in the table — not a dispatch layer.

| Keyword | Default (actual distribution) | Role |
|:--------|:------------------------------|:-----|
| `v_trim_5p` | `DiscreteUniform(0, 0)` | V 5′ trim (nt) |
| `v_trim_3p` | `DiscreteUniform(0, 10)` | V 3′ trim (nt) |
| `d_trim_5p` | `Geometric(0.20)` | D 5′ trim (nt; mean = 4) |
| `d_trim_3p` | `Geometric(0.21)` | D 3′ trim (nt; mean ≈ 3.8) |
| `j_trim_5p` | `DiscreteUniform(0, 14)` | J 5′ trim (nt) |
| `n1_length` | `DiscreteUniform(0, 18)` | N1 addition length (nt) |
| `n2_length` | `DiscreteUniform(0, 18)` | N2 addition length (nt) |
| `include_d` | `Bernoulli(0.99)` | rare VJ-only skip |
| `flank_5p` | `DomainFlankMix(...)` | 5′ unread flank |
| `flank_3p` | `DomainFlankMix(...)` | 3′ unread flank |
| `body_error_rate` | `GatedRate(0.40, Uniform(0.005, 0.07))` | SHM (`shm_igm`) |
| `indel_rate` | `GatedRate(0.08, Uniform(0.0, 0.002))` | SHM-like indels |
| `illumina_error` | `IlluminaError(0.001)` | instrument (`illumina_miseq`) |
| `min_length` | `80` | min accepted length (nt) |
| `max_length` | `900` | max accepted length (nt) |
| `max_retries` | `32` | recombine attempts before error |

```julia
using Distributions
params = train_params(;
    flank_5p = DomainFlankMix(DiscreteUniform(5, 40), DiscreteUniform(40, 120), 0.6),
    body_error_rate = shm_igg(),
    illumina_error = IlluminaError(0.001),  # merged PE Q30; alias illumina_miseq()
)
gen = ReadGenerator(db; params)
```

### Error model (SHM vs Illumina)

Both layers touch **VDJ body bases only** (V, N1, D, N2, J). Flanks are already
random DNA, so extra substitutions there would not change the training signal.
Order: recombination → SHM (+ rare indels) → Illumina substitutions.

**SHM** is zero-inflated, because a real IgM library is mostly germline with a
mutated memory tail. `GatedRate(p, rate)` does:

1. With probability `p`, draw a per-base substitution rate from `rate`.
2. Otherwise the rate is `0` (unmutated read).
3. Each body base then flips independently at that rate.

| Alias | Implementation | Meaning |
|:------|:---------------|:--------|
| `shm_none()` | `GatedRate(0.0, Dirac(0.0))` | naive, no SHM |
| `shm_igm()` | `GatedRate(0.40, Uniform(0.005, 0.07))` | ~60% germline, 40% at 0.5–7% |
| `shm_igg()` | `GatedRate(0.97, LogNormal(-2.65, 0.50))` | ~97% mutated, median ≈ 7%, right tail |
| `shm(p, lo, hi)` | `GatedRate(p, Uniform(lo, hi))` | custom |

This is a **training prior**, not a fit to one AIRR file. IgM uses a flat
Uniform over mutation depth; IgG uses `LogNormal` so high-SHM clones are not
clipped. We do **not** model AID hotspots (RGYW) or CDR vs framework targeting;
i.i.d. substitutions are a harder retrieval task (mutations everywhere).

**Illumina** is always-on and runs after SHM. Default `IlluminaError(0.001)` is
uniform ~Q30 (`0.1%`) on the merged contig — the right cartoon for **paired-end
MiSeq + overlap merge**. In the overlap, assemblers keep the higher-Q base, so
a single-mate error is usually dropped rather than randomly kept; the
single-end 3′ decay is cancelled. Substitutions only (Illumina indels ~10⁻⁶).
For unmerged single-end, use a ramp e.g. `IlluminaError(0.001, 0.004)`.
Disable with `IlluminaError(0.0)` / `illumina_none()`.

## License

MIT — see [LICENSE](LICENSE).

