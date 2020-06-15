using SFGAnalysis
using Documenter

makedocs(;
    modules=[SFGAnalysis],
    authors="Michael Lackner <lacknersmichael@gmail.com",
    repo="https://github.com/MLackner/SFGAnalysis.jl/blob/{commit}{path}#L{line}",
    sitename="SFGAnalysis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
deploydocs(
    repo   = "github.com/MLackner/SFGAnalysis.jl",
    branch = "gh-pages",
    devbranch = "master",
    devurl = "dev",
)