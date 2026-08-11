# germline.jl — Allele records and FASTA-backed germline DB.

struct Allele
    name::String
    sequence::String
end

Base.length(a::Allele) = length(a.sequence)
sequence_string(a::Allele) = a.sequence

function gene_from_allele(name::AbstractString)
    s = String(name)
    i = findfirst('*', s)
    isnothing(i) ? s : s[1:prevind(s, i)]
end

struct GermlineDB
    v::Vector{Allele}
    d::Vector{Allele}
    j::Vector{Allele}
end

n_v(db::GermlineDB) = length(db.v)
n_d(db::GermlineDB) = length(db.d)
n_j(db::GermlineDB) = length(db.j)

"""Subset alleles whose names are in `names` (exact match)."""
function filter_alleles(alleles::Vector{Allele}, names)
    want = names isa AbstractSet ? names : Set(String.(names))
    Allele[a for a in alleles if a.name in want]
end

"""Alleles whose names are **not** in `names`."""
function exclude_alleles(alleles::Vector{Allele}, names)
    drop = names isa AbstractSet ? names : Set(String.(names))
    Allele[a for a in alleles if !(a.name in drop)]
end

"""
    holdout_alleles(db, held_v, held_d, held_j) -> (train_db, held_db)

Split a germline for open-set protocols. Training streams use `train_db`
(held alleles removed). Evaluation galleries use `merge_germline(train, held)`
or the original full DB. Names are exact allele strings.
"""
function holdout_alleles(db::GermlineDB,
                         held_v = String[],
                         held_d = String[],
                         held_j = String[])
    hv = filter_alleles(db.v, held_v)
    hd = filter_alleles(db.d, held_d)
    hj = filter_alleles(db.j, held_j)
    train = GermlineDB(exclude_alleles(db.v, held_v),
                       exclude_alleles(db.d, held_d),
                       exclude_alleles(db.j, held_j))
    held = GermlineDB(hv, hd, hj)
    train, held
end

"""Concatenate allele lists (train ∪ held). Does not deduplicate names."""
merge_germline(a::GermlineDB, b::GermlineDB) =
    GermlineDB(vcat(a.v, b.v), vcat(a.d, b.d), vcat(a.j, b.j))

function read_fasta_alleles(path::AbstractString)
    alleles = Allele[]
    name = ""; buf = IOBuffer()
    open(path) do io
        for line in eachline(io)
            if startswith(line, '>')
                if !isempty(name)
                    push!(alleles, Allele(name, String(take!(buf))))
                end
                name = String(strip(line[2:end]))
                # drop description after whitespace
                sp = findfirst(isspace, name)
                !isnothing(sp) && (name = name[1:prevind(name, sp)])
            else
                print(buf, uppercase(strip(line)))
            end
        end
    end
    !isempty(name) && push!(alleles, Allele(name, String(take!(buf))))
    alleles
end

"""
    load_germline(v, d, j)
    load_germline(; v, d, j)

Build a [`GermlineDB`](@ref) from FASTA paths. `d` may be `nothing` / omitted
for VJ-only repertoires. Callers supply paths — no dataset naming is assumed.
"""
function load_germline(v_path::AbstractString, d_path::Union{AbstractString,Nothing},
                       j_path::AbstractString)
    v = read_fasta_alleles(v_path)
    j = read_fasta_alleles(j_path)
    d = isnothing(d_path) ? Allele[] : read_fasta_alleles(d_path)
    (isempty(v) || isempty(j)) && throw(ArgumentError("V and J FASTA must be non-empty"))
    GermlineDB(v, d, j)
end

load_germline(; v::AbstractString, j::AbstractString,
              d::Union{AbstractString,Nothing} = nothing) =
    load_germline(v, d, j)

sample_allele(rng::AbstractRNG, alleles::Vector{Allele}) =
    alleles[rand(rng, 1:length(alleles))]

function allele_gene_ids(alleles::Vector{Allele})
    genes = Dict{String,Int}()
    ids = Vector{Int}(undef, length(alleles))
    next = 1
    @inbounds for i in eachindex(alleles)
        g = gene_from_allele(alleles[i].name)
        if isempty(g)
            ids[i] = 0
        else
            ids[i] = get!(genes, g) do
                x = next
                next += 1
                x
            end
        end
    end
    ids
end

function allele_hamming_left(a::AbstractString, b::AbstractString)
    L = min(length(a), length(b))
    d = abs(length(a) - length(b))
    @inbounds for i in 1:L
        d += (a[i] != b[i])
    end
    d
end

"""Nearest same-gene sisters by left-aligned Hamming `(k, n)` (0 = pad)."""
function nearest_sister_indices(alleles::Vector{Allele}; k::Integer = 4)
    n = length(alleles)
    kk = max(Int(k), 0)
    out = zeros(Int, kk, n)
    kk == 0 && return out
    n < 2 && return out
    gene_ids = allele_gene_ids(alleles)
    by_gene = Dict{Int,Vector{Int}}()
    @inbounds for i in 1:n
        g = gene_ids[i]
        g == 0 && continue
        push!(get!(() -> Int[], by_gene, g), i)
    end
    for (_, members) in by_gene
        length(members) < 2 && continue
        for i in members
            dists = Tuple{Int,Int}[]
            ai = alleles[i].sequence
            for j in members
                j == i && continue
                push!(dists, (allele_hamming_left(ai, alleles[j].sequence), j))
            end
            sort!(dists; by = first)
            for r in 1:min(kk, length(dists))
                out[r, i] = dists[r][2]
            end
        end
    end
    out
end
