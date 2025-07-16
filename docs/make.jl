using Documenter, SymOrb

makedocs(format=Documenter.HTML(ansicolor=true,  assets = ["assets/style.css"]), sitename="SymOrb.jl", remotes=nothing,
pages = [
    "Home" => "index.md",
    "installation.md",
    "Basic usage" => [
        "usage/problem_definition.md",
        "usage/new_orbit.md",
        "usage/visualization.md"
        ],
    "Examples" => [
        "examples/examples.md",
        "examples/example_1.md",
        "examples/example_2.md",
        "examples/example_3.md",
        "examples/example_4.md",
        "examples/example_5.md",
        "examples/example_6.md",

    ],
    "APIs.md"
]
)

deploydocs(
    repo = "github.com/DipMathUniTO/SymOrb.jl.git",
)