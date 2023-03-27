push!(LOAD_PATH,"../src/")
using IFSintegrals
using Documenter

makedocs(
         sitename = "IFSintegrals.jl",
         modules  = [IFSintegrals],
         pages=[
                "Home" => "index.md",
                "Constructing fractals" => "makeIFS.md",
                "Plotting" => "plotting.md",
                "Approximating integrals" => "quadrature.md",
                "Integral equations" => "SIOs.md"
               ])
deploydocs(;
    repo="github.com/AndrewGibbs/IFSintegrals.jl",
)
