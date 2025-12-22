using PlasmaBO
using Documenter

DocMeta.setdocmeta!(PlasmaBO, :DocTestSetup, :(using PlasmaBO); recursive = true)

makedocs(;
    modules = [PlasmaBO],
    authors = "Beforerr <zzj956959688@gmail.com> and contributors",
    sitename = "PlasmaBO.jl",
    format = Documenter.HTML(;
        canonical = "https://JuliaSpacePhysics.github.io/PlasmaBO.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md"
    ],
)

deploydocs(;
    repo = "github.com/JuliaSpacePhysics/PlasmaBO.jl",
    push_preview = true,
)
