using Documenter
using UniqueKronecker
using DocumenterCitations

ENV["JULIA_DEBUG"] = "Documenter"
DocMeta.setdocmeta!(UniqueKronecker, :DocTestSetup, :(using UniqueKronecker); recursive=true)

PAGES = [
    "Home" => "index.md",
    "Unique Kronecker Products" => [
        "Basics" => "uniquekronecker/basics.md",
        "Higher Order" => "uniquekronecker/higherorder.md",
    ],
    "Circulant Kronecker Products" => "circulantkronecker/basics.md",
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
    sitename = "Unique Kronecker",
    clean = true, doctest = false, linkcheck = false,
    authors = "Tomoki Koike <tkoike45@gmail.com>",
    repo = Remotes.GitHub("smallpondtom", "UniqueKronecker.jl"),
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = "https://github.com/smallpondtom/UniqueKronecker.jl",
        assets=String[
            "assets/citations.css",
            "assets/favicon.ico",
        ],
        # analytics = "G-B2FEJZ9J99",
    ),
    modules = [UniqueKronecker,],
    pages = PAGES,
    plugins=[bib],
)

deploydocs(
    repo = "github.com/smallpondtom/UniqueKronecker.jl",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
    # Add other deployment options as needed
)