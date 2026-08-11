# regions.jl — Region labels and spans.

@enum RegionLabel::UInt8 begin
    REG_PAD = 0
    REG_V = 1
    REG_N1 = 2
    REG_D = 3
    REG_N2 = 4
    REG_J = 5
end

struct Span
    start::Int
    stop::Int
end

const EMPTY_SPAN = Span(0, -1)
Base.isempty(s::Span) = s.stop < s.start
Base.length(s::Span) = isempty(s) ? 0 : s.stop - s.start + 1

"""Longest contiguous run of each locus label → named spans."""
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
    (flank5 = EMPTY_SPAN,  # flanks are PAD; use geometry from assemble metadata
     v = pick(REG_V),
     n1 = pick(REG_N1),
     d = pick(REG_D),
     n2 = pick(REG_N2),
     j = pick(REG_J),
     flank3 = EMPTY_SPAN)
end
