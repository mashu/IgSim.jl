"""
    IgSim

Streaming V(D)J read simulator for training. Uniform allele sampling; defaults
via [`train_params`](@ref). No gold-identity gate.
"""
module IgSim

using Random
using Distributions
using Statistics
using CodecZlib

include("alphabet.jl")
include("regions.jl")
include("germline.jl")
include("params.jl")
include("trim.jl")
include("nadd.jl")
include("recombine.jl")
include("assemble.jl")
include("noise.jl")
include("read.jl")
include("stream.jl")
include("contrastive.jl")
include("batch.jl")

export
    NucleotideAlphabet, tokenize, detokenize, DNA_BASES,
    RegionLabel, REG_PAD, REG_V, REG_N1, REG_D, REG_N2, REG_J,
    Span, EMPTY_SPAN, spans_from_labels,
    Allele, GermlineDB, load_germline, read_fasta_alleles,
    n_v, n_d, n_j, sequence_string,
    gene_from_allele, allele_gene_ids, nearest_sister_indices,
    holdout_alleles, merge_germline, filter_alleles, exclude_alleles,
    HoldoutVGenerator, MixedReadGenerator, mixed_curriculum,
    SimParams, DomainFlankMix, GatedRate, train_params,
    default_flank_5p, default_flank_3p,
    easy_params, mid_params, hard_params, short_flank_params,
    LabeledRead, ReadGenerator,
    ReadStream, take_batch, channel_stream,
    ContrastiveExample, ContrastiveStream,
    collate_reads, ReadBatch

end # module
