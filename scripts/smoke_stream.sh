#!/usr/bin/env bash
# Smoke against test fixtures (no external datasets).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
julia --project=. -e '
using IgSim, Random, Statistics
fix = joinpath(@__DIR__, "test", "fixtures")
db = load_germline(; v=joinpath(fix,"V.fasta"), d=joinpath(fix,"D.fasta"), j=joinpath(fix,"J.fasta"))
println("germline V=$(n_v(db)) D=$(n_d(db)) J=$(n_j(db))")
gen = ReadGenerator(db)
batch = take_batch(gen, 128; seed=7)
println("batch=$(length(batch)) mean_len=$(round(mean(length(r.sequence) for r in batch); digits=1))")
println("mean_flank5=$(round(mean(r.flank5 for r in batch); digits=1)) mean_flank3=$(round(mean(r.flank3 for r in batch); digits=1))")
cs = ContrastiveStream(gen; sister_k=2, seed=9)
cex = take_batch(cs, 4)
for (i, e) in enumerate(cex)
    println("contrastive[$i] V=$(e.read.v_call) pos=$(e.positive_v) n_neg=$(length(e.negative_v))")
end
t0 = time(); n = 2000; take_batch(gen, n; seed=1); dt = time() - t0
println("throughput ≈ $(round(n/dt; digits=0)) reads/s")
'
