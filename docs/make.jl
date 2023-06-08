push!(LOAD_PATH,"../src/")
using IFSintegrals
using Documenter

# make logos
include("make_logo.jl")
# for light backgrounds
make_logo("logo";line_colour="black")
# for dark backgrounds
make_logo("logo-dark";line_colour="white")

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
