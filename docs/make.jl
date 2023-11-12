push!(LOAD_PATH,"../src/")
using IFSintegrals
using Documenter

makedocs(
         sitename = "IFSintegrals.jl",
         modules  = [IFSintegrals],
         checkdocs=:none,
         pages=[
                "Home" => "index.md",
                "Constructing fractals" => "makeIFS.md",
                "Plotting" => "plotting.md",
                "Approximating integrals" => "quadrature.md",
                "Integral equations" => "SIOs.md"
               ],
               format = Documenter.HTML(
                prettyurls = get(ENV, "CI", nothing) == "true"))
deploydocs(;
    repo="github.com/AndrewGibbs/IFSintegrals.jl",
)
