push!(LOAD_PATH,"../src/")
using IFSintegrals
using Documenter

makedocs(
         sitename = "IFSintegrals.jl",
         modules  = [IFSintegrals],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/AndrewGibbs/IFSintegrals.jl",
)
