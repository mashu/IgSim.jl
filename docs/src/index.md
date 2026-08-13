# IgSim.jl

Streaming V(D)J read simulator for training: uniform allele sampling, labeled
spans, holdout splits, and contrastive sister negatives.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/mashu/IgSim.jl")
```

## Quick start

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
train_gen = ReadGenerator(train_db)                  # never samples held V
held_gen  = HoldoutVGenerator(full_db, held_v_names) # gold V ∈ held set
# evaluation gallery = full_db (or merge_germline(train_db, held))
```

### Custom parameters

Defaults come from [`train_params`](@ref). Every keyword is overridable — see
[Parameters](@ref) for the full table and defaults.

```julia
using Distributions
params = train_params(;
    flank_5p = DomainFlankMix(DiscreteUniform(5, 40), DiscreteUniform(40, 120), 0.6),
    body_error_rate = GatedRate(0.97, Uniform(0.015, 0.13)),  # IgG; alias shm_igg()
    illumina_error = IlluminaError(0.001),
    include_d = Bernoulli(1.0),
)
gen = ReadGenerator(db; params)
```

