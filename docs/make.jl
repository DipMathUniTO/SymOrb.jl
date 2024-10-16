using Documenter, SymOrb

makedocs(format=Documenter.HTML(ansicolor=true), sitename="SymOrb.jl", remotes=nothing,
pages = [
    "index.md",
    "installation.md",
    "Basic usage" => [
        "usage/problem_definition.md",
        "usage/new_orbit.md",
        "usage/visualization.md"
        ],
    "APIs.md"
]
)

deploydocs(
    repo = "github.com/DipMathUniTO/SymOrb.jl.git",
)