__precompile__()

include("src/IFSintegrals.jl")
using .IFSintegrals
using MAT, Printf, FastGaussQuadrature
using FileIO, JLD2
## ------------------------------------ set the parameters -------------------------------- ##
# K = [3, 23, 53]
# A = [1/3, 0.2]
α = 0.0#parse(Float64,ARGS[2])
ℓ_min = 6#parse(Int64,ARGS[3])
ℓ_max = 6#parse(Int64,ARGS[4])
# k = 2.0
# α = 1/3
# ℓ_min = 1
# ℓ_max = 3
C_guess = 1.0000001
# d = [0 1/sqrt(2) -1/sqrt(2)]

## ---- First define a load of functions which will get used in the convergence test -------- ##

# get_quad_from_fda(C::Float64,h_bem::Float64,d::Float64) = C*h_bem^(1+d/2)
get_quad_from_fda(C::Float64,ℓ::Int64,Γ::Attractor) = C_guess*Γ.diameter*Γ.IFS[1].r^(ℓ_max)#C*(Γ.diameter*Γ.IFS[1].r^ℓ)^(1+Γ.Hausdorff_dimension/2)

function CantorDustConvergence(α::Float64, C_guess::Float64, ℓ::Int64)
    Γ = CantorDust(α)
    Cosc = Inf
    ρ = Γ.IFS[1].r
    h_BEM = 1.0000000000001*Γ.diameter*ρ^ℓ
    
    # solve the scattering problem, get the solution on the fractal
    # h_quad = get_quad_from_fda(C_guess,h_BEM,Γ.Hausdorff_dimension)
    h_quad = get_quad_from_fda(C_guess,ℓ,Γ)

    Si = SingleLayer(Γ, im)
    Sih = DiscreteBIO(Si,h_BEM,h_quad)
    Gi = Sih.Galerkin_matrix

    ℓ = length(Sih.Lₕ[1])
    FileIO.save(@sprintf("sols/sol_ki_ell%d_gap_%.2f.jld2",ℓ,α),
    "gap",α,"C_guess",C_guess,"h_BEM",h_BEM,"h_quad",h_quad,
    "ell",ℓ,"Gi",Gi,"Cosc",Cosc,"Lh",Sih.Lₕ,"k",im)
end

## ----------------------------- Now run for the different parameters -------------------- ##

for ℓ =ℓ_min:ℓ_max
    CantorDustConvergence(α, C_guess, ℓ)
end

# https://discourse.julialang.org/t/save-variable/33986
