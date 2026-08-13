# batch.jl — Pad collate for downstream training (no Flux/device).

"""Padded batch of [`LabeledRead`](@ref) (sequences, region labels, allele calls)."""
struct ReadBatch
    sequences::Vector{String}
    labels::Matrix{UInt8}          # (T, B) padded with REG_PAD
    lengths::Vector{Int32}
    v_call::Vector{String}
    d_call::Vector{Union{String,Missing}}
    j_call::Vector{String}
    multi_d::Vector{Bool}
    reads::Vector{LabeledRead}
end

"""Pad-collate a vector of reads into a [`ReadBatch`](@ref)."""
function collate_reads(reads::Vector{LabeledRead})
    B = length(reads)
    B == 0 && return ReadBatch(String[], Matrix{UInt8}(undef, 0, 0), Int32[],
                               String[], Union{String,Missing}[], String[],
                               Bool[], LabeledRead[])
    lens = Int32[length(r.sequence) for r in reads]
    T = Int(maximum(lens))
    labels = fill(UInt8(REG_PAD), T, B)
    @inbounds for b in 1:B
        L = Int(lens[b])
        labels[1:L, b] .= reads[b].labels
    end
    ReadBatch([r.sequence for r in reads], labels, lens,
              [r.v_call for r in reads],
              Union{String,Missing}[r.d_call for r in reads],
              [r.j_call for r in reads],
              [r.multi_d for r in reads],
              reads)
end
