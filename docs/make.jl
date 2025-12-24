using PlasmaBO
using Documenter, DocumenterCitations

DocMeta.setdocmeta!(PlasmaBO, :DocTestSetup, :(using PlasmaBO); recursive = true)

const bib = CitationBibliography(joinpath(@__DIR__, "PlasmaBO.jl.bib"))

makedocs(;
    modules = [PlasmaBO],
    authors = "Beforerr <zzj956959688@gmail.com> and contributors",
    sitename = "PlasmaBO.jl",
    format = Documenter.HTML(;
        canonical = "https://JuliaSpacePhysics.github.io/PlasmaBO.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Firehose instability" => "firehose_Astfalk17.md",
        "Ion cyclotron emission" => "ice_Irvine18.md",
    ],
    plugins = [bib],
)

deploydocs(;
    repo = "github.com/JuliaSpacePhysics/PlasmaBO.jl",
    push_preview = true,
)
