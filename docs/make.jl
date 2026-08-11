using Documenter
using IgSim

DocMeta.setdocmeta!(IgSim, :DocTestSetup, :(using IgSim); recursive = true)

makedocs(;
    modules = [IgSim],
    authors = "Mateusz Kaduk",
    sitename = "IgSim.jl",
    repo = Remotes.GitHub("mashu", "IgSim.jl"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://mashu.github.io/IgSim.jl",
        edit_link = "main",
    ),
    pages = [
        "Home" => "index.md",
        "Parameters" => "parameters.md",
        "API" => "api.md",
    ],
    checkdocs = :exports,
    warnonly = [:missing_docs],
)

deploydocs(;
    repo = "github.com/mashu/IgSim.jl",
    devbranch = "main",
    push_preview = true,
)
