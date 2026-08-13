using Test
using Random
using Distributions
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
    @test train_params().body_error_rate isa GatedRate
end

@testset "gated SHM + curriculum" begin
    rng = MersenneTwister(7)
    @test all(==(0.0), [rand(rng, GatedRate(0.0, Dirac(0.1))) for _ in 1:40])
    @test all(==(0.07), [rand(rng, GatedRate(1.0, Dirac(0.07))) for _ in 1:20])
    @test shm_none().p == 0
    @test shm_igm().p == 0.40
    @test shm_igg().p == 0.97
    @test shm_igg().rate isa LogNormal
    @test train_params().body_error_rate.p == shm_igm().p
    @test hard_params().body_error_rate.p == shm_igg().p
    @test all(==(0.04), [rand(rng, shm(1.0, 0.04, 0.04)) for _ in 1:20])
    db = load_germline(; v = joinpath(FIX, "V.fasta"),
                         d = joinpath(FIX, "D.fasta"),
                         j = joinpath(FIX, "J.fasta"))
    easy = ReadGenerator(db; params = easy_params())
    @test all(r -> r.body_error_rate == 0 && r.n_errors == 0,
              [easy(rng) for _ in 1:30])
    mix = mixed_curriculum(db)
    @test mix isa MixedReadGenerator
    r = mix(rng)
    @test !isempty(r.sequence)
    @test length(r.sequence) == length(r.labels)
end

@testset "D remnant keeps label; empty D is no-D" begin
    db = load_germline(; v = joinpath(FIX, "V.fasta"),
                         d = joinpath(FIX, "D.fasta"),
                         j = joinpath(FIX, "J.fasta"))
    d1 = filter_alleles(db.d, ["IGHD1-1*01"])
    @test length(d1) == 1
    @test length(d1[1].sequence) == 25
    db1 = GermlineDB(db.v, d1, db.j)
    rng = MersenneTwister(11)
    keep = train_params(include_d = Bernoulli(1.0),
                        d_trim_5p = DiscreteUniform(12, 12),
                        d_trim_3p = DiscreteUniform(12, 12))
    parts = IgSim.recombine(rng, db1, keep)
    @test !isnothing(parts)
    @test parts.has_d
    @test parts.d_call == "IGHD1-1*01"
    @test length(parts.d_body) == 1
    eaten = train_params(include_d = Bernoulli(1.0),
                         d_trim_5p = DiscreteUniform(20, 20),
                         d_trim_3p = DiscreteUniform(20, 20))
    parts0 = IgSim.recombine(MersenneTwister(11), db1, eaten)
    @test !isnothing(parts0)
    @test !parts0.has_d
    @test ismissing(parts0.d_call)
    @test isempty(parts0.d_body)
    @test train_params().d_trim_5p isa Geometric
    @test train_params().d_trim_3p isa Geometric
    rng2 = MersenneTwister(21)
    gen = ReadGenerator(db; params = train_params())
    n_miss = count(r -> !r.has_d, [gen(rng2) for _ in 1:250])
    @test n_miss < 80
end

@testset "tandem DD: DJ then D+DJ then V" begin
    db = load_germline(; v = joinpath(FIX, "V.fasta"),
                         d = joinpath(FIX, "D.fasta"),
                         j = joinpath(FIX, "J.fasta"))
    @test n_d(db) == 2
    @test train_params().include_dd isa Bernoulli
    @test train_params().include_dd.p == 0.001
    off = train_params(include_d = Bernoulli(1.0),
                       include_dd = Bernoulli(0.0))
    @test all(r -> !r.multi_d, [ReadGenerator(db; params = off)(MersenneTwister(i))
                                for i in 1:40])
    on = train_params(include_d = Bernoulli(1.0),
                      include_dd = Bernoulli(1.0),
                      d_trim_5p = DiscreteUniform(2, 2),
                      d_trim_3p = DiscreteUniform(2, 2),
                      j_trim_5p = DiscreteUniform(0, 0),
                      n2_length = DiscreteUniform(3, 3),
                      body_error_rate = shm_none(),
                      indel_rate = GatedRate(0.0, Dirac(0.0)),
                      illumina_error = illumina_none())
    rng = MersenneTwister(13)
    parts = IgSim.recombine(rng, db, on)
    @test !isnothing(parts)
    @test parts.has_d
    @test parts.multi_d
    @test parts.d_call isa String
    names = split(parts.d_call, ',')
    @test length(names) == 2
    @test names[1] != names[2]
    @test Set(names) == Set(a.name for a in db.d)
    @test !occursin(' ', parts.d_call)
    @test length(parts.d_body) == length(db.d[findfirst(a -> a.name == names[1], db.d)].sequence) - 4
    @test length(parts.d2_body) == length(db.d[findfirst(a -> a.name == names[2], db.d)].sequence) - 4
    @test length(parts.n3) == 3
    @test length(parts.n2) == 3
    gen = ReadGenerator(db; params = on)
    r = gen(MersenneTwister(13))
    @test r.multi_d === true
    @test r.has_d
    @test occursin(',', r.d_call)
    @test any(==(UInt8(REG_N3)), r.labels)
    @test !isempty(r.spans.d)
    @test !isempty(r.spans.d2)
    @test !isempty(r.spans.n3)
    @test r.spans.d.stop < r.spans.n3.start
    @test r.spans.n3.stop < r.spans.d2.start
    rb = collate_reads([r])
    @test rb.multi_d == [true]
    eaten5 = train_params(include_d = Bernoulli(1.0),
                          include_dd = Bernoulli(1.0),
                          d_trim_5p = DiscreteUniform(40, 40),
                          d_trim_3p = DiscreteUniform(0, 0),
                          n2_length = DiscreteUniform(2, 2))
    # First DJ D is fully eaten → no D, so no DD either.
    parts0 = IgSim.recombine(MersenneTwister(5), db, eaten5)
    @test !isnothing(parts0)
    @test !parts0.has_d
    @test !parts0.multi_d
    @test ismissing(parts0.d_call)
    one_d = GermlineDB(db.v, filter_alleles(db.d, ["IGHD1-1*01"]), db.j)
    parts1 = IgSim.recombine(MersenneTwister(5), one_d, on)
    @test !isnothing(parts1)
    @test parts1.has_d
    @test !parts1.multi_d
    @test parts1.d_call == "IGHD1-1*01"
    # Short D eaten, long D may remain: never gold DD.
    mid_eat = train_params(include_d = Bernoulli(1.0),
                           include_dd = Bernoulli(1.0),
                           d_trim_5p = DiscreteUniform(26, 26),
                           d_trim_3p = DiscreteUniform(0, 0),
                           j_trim_5p = DiscreteUniform(0, 0))
    for seed in 1:30
        p = IgSim.recombine(MersenneTwister(seed), db, mid_eat)
        isnothing(p) && continue
        @test !p.multi_d
        if p.has_d
            @test p.d_call == "IGHD2-2*01"
        end
    end
end

@testset "trim/N helpers" begin
    @test IgSim.trim_5p_ascii("ACGT", 1) == "CGT"
    @test IgSim.trim_3p_ascii("ACGT", 1) == "ACG"
    @test IgSim.sample_n(MersenneTwister(0), 0) == ""
    @test length(IgSim.sample_n(MersenneTwister(0), 5)) == 5
end

@testset "Illumina body-only substitutions" begin
    @test train_params().illumina_error == illumina_miseq()
    @test illumina_miseq().p5 == illumina_miseq().p3 == 0.001
    @test illumina_none().p5 == 0 && illumina_none().p3 == 0
    seq = "AAAAACGTAAAA"
    labs = vcat(fill(UInt8(REG_PAD), 4),
                fill(UInt8(REG_V), 4),
                fill(UInt8(REG_PAD), 4))
    rng = MersenneTwister(3)
    out, n = IgSim.apply_illumina_error(rng, seq, labs, IlluminaError(1.0))
    @test out[1:4] == "AAAA"
    @test out[9:12] == "AAAA"
    @test out[5:8] != "ACGT"
    @test n == 4
    none, n0 = IgSim.apply_illumina_error(MersenneTwister(3), seq, labs,
                                          illumina_none())
    @test none == seq
    @test n0 == 0
    db = load_germline(; v = joinpath(FIX, "V.fasta"),
                         d = joinpath(FIX, "D.fasta"),
                         j = joinpath(FIX, "J.fasta"))
    gen = ReadGenerator(db; params = train_params(
        body_error_rate = shm_none(),
        indel_rate = GatedRate(0.0, Dirac(0.0)),
        illumina_error = illumina_none()))
    r = gen(MersenneTwister(1))
    @test r.n_errors == 0
    @test r.n_illumina == 0
end
