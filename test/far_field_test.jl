using IFSintegrals, MAT, LinearAlgebra, Test

α = 1/3
Γ = CantorSet(contraction=α);
k = 30.0
d = [0.5000, -0.8660]

ℓ = 10
thresh = 1E-3

f(x::Float64) = exp(im*k*d[1]*x)
Sₖ = SingleLayer(Γ,k)

h_θ = 2π/300
θ = 0:h_θ:2π
FFvecs = zeros(Complex{Float64},length(θ))

h_BEM = α^ℓ*(1+10*eps())
h_quad = h_BEM*α^2
Sₖₕ = DiscreteSIO(Sₖ; h_mesh = h_BEM, h_quad = h_quad);
ϕₕ = Sₖₕ\f
Fϕₕ = far_field_pattern(ϕₕ,k; h_quad = h_quad)
FF_const = -sqrt(1im/(8*π*k))
FF_Hausdorff = FF_const*Fϕₕ.(θ)

# AM = matread("ff_ref.mat")
# FF_Lebesgue = AM["ExperimentVec"][14]["FarField"]

AM2 = matread("Lebesgue_FF_vals.mat")
FF_Lebesgue = AM2["FF_vals"]

FF_Lebesgue = reshape(FF_Lebesgue,length(FF_Lebesgue))
FF_err = norm(FF_Hausdorff-FF_Lebesgue)/norm(FF_Lebesgue)

@test(FF_err<thresh)