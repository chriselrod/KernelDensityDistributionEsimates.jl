using Documenter, KernelDensityDistributionEsimates

makedocs(;
    modules=[KernelDensityDistributionEsimates],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/chriselrod/KernelDensityDistributionEsimates.jl/blob/{commit}{path}#L{line}",
    sitename="KernelDensityDistributionEsimates.jl",
    authors="Chris Elrod",
    assets=[],
)

deploydocs(;
    repo="github.com/chriselrod/KernelDensityDistributionEsimates.jl",
)
