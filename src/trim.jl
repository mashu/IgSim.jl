# trim.jl — Exonuclease trims (ASCII DNA alleles).

trim_5p_ascii(seq::String, n::Integer) =
    (n = max(Int(n), 0); n >= length(seq) ? "" : seq[n + 1:end])

trim_3p_ascii(seq::String, n::Integer) =
    (n = max(Int(n), 0); n >= length(seq) ? "" : seq[1:length(seq) - n])
