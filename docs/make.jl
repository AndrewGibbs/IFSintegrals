push!(LOAD_PATH,"../src/")
using IFSintegrals Documentermakedocs(
         sitename = "IFSintegrals.jl",
         modules  = [IFSintegrals],
         pages=[
                "Home" => "index.md"
               ])deploydocs(;
    repo="github.com/AndrewGibbs/IFSintegrals.jl",
)
