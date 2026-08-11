using Test
using Random
using IgSim

const FIX = joinpath(@__DIR__, "fixtures")

@testset "IgSim smoke" begin
    db = load_germline(; v = joinpath(FIX, "V.fasta"),
                         d = joinpath(FIX, "D.fasta"),
                         j = joinpath(FIX, "J.fasta"))
    @test n_v(db) == 3
    @test n_j(db) == 2
    gen = ReadGenerator(db; params = train_params())
    rng = MersenneTwister(1)
    r = gen(rng)
    @test !isempty(r.sequence)
    @test length(r.sequence) == length(r.labels)
    @test !isempty(r.v_call)
    @test !isempty(r.j_call)
    @test !isempty(r.spans.v)
    batch = take_batch(gen, 32; seed = 2)
    @test length(batch) == 32
    rb = collate_reads(batch)
    @test size(rb.labels, 2) == 32
    cs = ContrastiveStream(gen; sister_k = 2, seed = 3)
    ce = first(iterate(cs))
    @test ce.positive_v > 0
    @test !(ce.positive_v in ce.negative_v)
end

@testset "holdout + domain flanks" begin
    db = load_germline(; v = joinpath(FIX, "V.fasta"),
                         d = joinpath(FIX, "D.fasta"),
                         j = joinpath(FIX, "J.fasta"))
    held_names = [db.v[1].name]
    train, held = holdout_alleles(db, held_names, String[], String[])
    @test n_v(train) == n_v(db) - 1
    @test n_v(held) == 1
    hg = HoldoutVGenerator(db, held_names)
    r = hg(MersenneTwister(1))
    @test r.v_call in Set(held_names)
    @test train_params().flank_5p isa DomainFlankMix
end

@testset "trim/N helpers" begin
    @test IgSim.trim_5p_ascii("ACGT", 1) == "CGT"
    @test IgSim.trim_3p_ascii("ACGT", 1) == "ACG"
    @test IgSim.sample_n(MersenneTwister(0), 0) == ""
    @test length(IgSim.sample_n(MersenneTwister(0), 5)) == 5
end
