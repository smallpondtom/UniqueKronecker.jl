using Documenter
using UniqueKronecker
using DocumenterCitations

ENV["JULIA_DEBUG"] = "Documenter"

PAGES = [
    "Home" => "index.md",
    "Unique Kronecker Products" => [
        "Basics" => "uniquekronecker/basics.md",
        "Higher Order" => "uniquekronecker/higherorder.md",
    ],
    "Matrix Conversions" => [
        "Polynomial Dynamical System" => "matrices/dynamicalsystem.md",
        "Elimination Matrix" => "matrices/elimination.md",
        "Duplication Matrix" => "matrices/duplication.md",
        "Commutation Matrix" => "matrices/commutation.md",
        "Symmetrizer Matrix" => "matrices/symmetrizer.md"
    ],
    "API Reference" => "api.md",
    "Paper Reference" => "paper.md",
]

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(
    sitename = "LiftAndLearn.jl",
    clean = true, doctest = false, linkcheck = false,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = "https://github.com/smallpondtom/UniqueKronecker.jl",
        assets=String[
            "assets/citations.css",
        ],
        # analytics = "G-B2FEJZ9J99",
    ),
    modules = [
        UniqueKronecker,
    ],
    pages = PAGES,
    plugins=[bib],
)

deploydocs(
    repo = "github.com/smallpondtom/UniqueKronecker.jl.git",
    branch = "gh-pages",
    devbranch = "main",
    # Add other deployment options as needed
)