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
        size_threshold_warn = 1024_000,
        size_threshold = 4096_000,
    ),
    pages = [
        "Home" => "index.md",
        "Ring beam instability" => "ringbeam_Umeda12.md",
        "Ion beam instability" => "ionbeam_gary84.md",
        "Firehose instability" => "firehose_Astfalk17.md",
        "Mirror mode" => "mirror_mode.md",
        "Ion cyclotron emission" => "ice_Irvine18.md",
        "R/L/P modes (PBK solver)" => "rlp_Cattaert07.md",
        "Cold plasma (fluid vs kinetic solver)" => "cold_plasma.md",
        "Dispersion surface tracking" => "dispersion_surface.md",
        "Wave polarization and handedness" => "demo_polarization.md",
        "Math notes" => "math_notes.md",
    ],
    plugins = [bib],
)

deploydocs(;
    repo = "github.com/JuliaSpacePhysics/PlasmaBO.jl",
    push_preview = true,
)
