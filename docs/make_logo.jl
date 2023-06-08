using IFSintegrals, Plots
import IFSintegrals: mesh_fractal

function make_logo(file_name;line_colour="black")
Γ = KochFlake()
h_mesh = 0.45#*sqrt(3)
𝙈 = mesh_fractal(Γ,h_mesh)
x = rand(length(𝙈))
xautolimitmargin = (0.0, 0.0)
yautolimitmargin = (0.0, 0.0)
p = plot(𝙈,x,levels=4,c=:Dark2_7,
    background=nothing,linecolor=line_colour,
    colorbar=nothing,axis=([], false),
    linewidth = 1.5, aspect_ratio=1.0,
    size = (875, 1000),
    ylim=(-1,1),
    xlim=(-0.875,0.875))
# display(p)
savefig(p,"src/assets/"*file_name*".png")
#:buda10, :viridis, Dark2_7,PRGn
end

# include("make_logo.jl")
# for light backgrounds
make_logo("logo";line_colour="black")
# for dark backgrounds
make_logo("logo-dark";line_colour="white")