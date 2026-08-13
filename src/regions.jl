# regions.jl — Region labels and spans.

@enum RegionLabel::UInt8 begin
    REG_PAD = 0
    REG_V = 1
    REG_N1 = 2
    REG_D = 3
    REG_N2 = 4
    REG_J = 5
    REG_N3 = 6  # N-addition between two Ds (tandem DD)
end

struct Span
    start::Int
    stop::Int
end

const EMPTY_SPAN = Span(0, -1)
Base.isempty(s::Span) = s.stop < s.start
Base.length(s::Span) = isempty(s) ? 0 : s.stop - s.start + 1

"""Contiguous runs of `lab` in 5′→3′ order."""
function ordered_runs(labels::AbstractVector{<:Integer}, lab)
    u = UInt8(lab)
    spans = Span[]
    i = 1
    n = length(labels)
    while i <= n
        if UInt8(labels[i]) == u
            j = i
            while j <= n && UInt8(labels[j]) == u
                j += 1
            end
            push!(spans, Span(i, j - 1))
            i = j
        else
            i += 1
        end
    end
    spans
end

"""Longest contiguous run of each locus label → named spans.

`d` / `d2` are the first and second D runs in 5′→3′ order (tandem DD).
Other loci still use the longest run.
"""
function spans_from_labels(labels::AbstractVector{<:Integer})
    best = Dict{UInt8,Tuple{Int,Int,Int}}()  # lab => (len, start, stop)
    i = 1
    n = length(labels)
    while i <= n
        lab = UInt8(labels[i])
        j = i
        while j <= n && UInt8(labels[j]) == lab
            j += 1
        end
        len = j - i
        prev = get(best, lab, (0, 0, -1))
        if len > prev[1]
            best[lab] = (len, i, j - 1)
        end
        i = j
    end
    function pick(lab)
        t = get(best, UInt8(lab), nothing)
        isnothing(t) && return EMPTY_SPAN
        Span(t[2], t[3])
    end
    ds = ordered_runs(labels, REG_D)
    d = isempty(ds) ? EMPTY_SPAN : ds[1]
    d2 = length(ds) >= 2 ? ds[2] : EMPTY_SPAN
    (flank5 = EMPTY_SPAN,  # flanks are PAD; use geometry from assemble metadata
     v = pick(REG_V),
     n1 = pick(REG_N1),
     d = d,
     n3 = pick(REG_N3),
     d2 = d2,
     n2 = pick(REG_N2),
     j = pick(REG_J),
     flank3 = EMPTY_SPAN)
end
