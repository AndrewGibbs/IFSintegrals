# using IFSintegrals, LinearAlgebra

# k = 10.0

function get_nf_and_ff(Γ,f,k; amb_dim = Γ.spatial_dimension)
    S = SingleLayerOperatorHelmholtz(Γ, k, ambient_dimension=amb_dim)
    Sₕ = DiscreteSIO(S, h_mesh = 0.1, h_quad = 0.01)
    ϕₕ = Sₕ \ f
    𝙎ϕₕ = SingleLayerPotentialHelmholtz(ϕₕ, k, ambient_dimension=amb_dim)
    uˢ = x->-1*𝙎ϕₕ(x)
    FFP = far_field_pattern(ϕₕ, k, ambient_dimension=amb_dim)
    return uˢ, FFP
end

# d = [0,-1,0];

# γ = CantorSet()
# Γ = CantorDust()

function ff_nf_agreement_R2_screen(Γ,r,d,k)
    f₁(x) = exp(im*k*(d[1]*x[1]))
    NF,FF = get_nf_and_ff(Γ, f₁, k; amb_dim=2)
    θ = rand()
    x = r*[cos(θ),sin(θ)]
    SH_1 = exp(im*k*r) * FF(θ) /sqrt(r)
    SH_2 = NF(x)
    return isapprox(SH_1,SH_2,rtol=10/r)
end

function ff_nf_agreement_R2(Γ,r,d,k)
    f₂(x) = exp(im*k*(d[1]*x[1]+d[2]*x[2]))
    NF,FF = get_nf_and_ff(Γ, f₂, k)
    θ = rand()
    x = r*[cos(θ),sin(θ)]
    SH_1 = exp(im*k*r) * FF(θ) /sqrt(r)
    SH_2 = NF(x)
    return isapprox(SH_1,SH_2,rtol=10/r)
end

function ff_nf_agreement_R3_screen(Γ,r,d,k)
    f₂(x) = exp(im*k*(d[1]*x[1]+d[2]*x[2]))
    NF,FF = get_nf_and_ff(Γ, f₂, k; amb_dim=3)
    θ = rand()
    φ = rand()
    x = r*[sin(θ)*cos(φ),sin(θ)*sin(φ),cos(θ)]
    SH_1 = exp(im*k*r) * FF([θ, φ]) /r
    SH_2 = NF(x)
    return isapprox(SH_1,SH_2,rtol=100/r)
end

# ff_nf_agreement_R2_screen(γ,1000,[0,-1],k)

# ff_nf_agreement_R2(Γ,1000,[0,-1],k)

# ff_nf_agreement_R3_screen(Γ,100000,[0,-1,0],k)