using Documenter, SymPath

makedocs(format=Documenter.HTML(ansicolor=true), sitename="SymPath.jl",
pages = [
    "index.md",
    "installation.md",
    "APIs.md"
    ]
)
