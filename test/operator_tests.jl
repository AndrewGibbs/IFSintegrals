# using IFSintegrals, LinearAlgebra, Test

function surface_test(Γ::SelfSimilarFractal)
    k = 1 + rand()*10
    f(x) = sin(norm(x)^2) # Laplace data

    amb_dim = Γ.spatial_dimension
    x = Vector{Float64}(Γ.barycentre) + Γ.diameter*(ones(amb_dim)+rand(amb_dim)) # random test point
    d = rand(amb_dim) # random 
    d ./ norm(d) # incident direction, normalised
    if amb_dim == 2
        θ = rand()
    elseif amb_dim == 2
        θ = rand(2)
    end

    g(x) = exp(im*k*(d'*x)) # pw data

    Sₖ = SingleLayerOperatorHelmholtz(Γ, k)
    Sₖₕ = DiscreteSIO(Sₖ; h_mesh=0.1, h_quad=0.05)
    ϕₖₕ = Sₖₕ\g
    𝙎ₖϕₕ = SingleLayerPotentialHelmholtz(ϕₖₕ, k; h_quad = 0.05)
    𝙎ₖϕₕ(x)

    Sₖ = SingleLayerOperatorHelmholtz(Γ, k)
    Sₖₕ = DiscreteSIO(Sₖ; h_mesh=0.1, h_quad=0.05)
    ϕₖₕ = Sₖₕ\g
    # 𝙁ₖϕₕ = far_field_pattern(ϕₖₕ,k; h_quad = 0.05)
    # 𝙁ₖϕₕ(θ)

    S₀ = SingleLayerOperatorLaplace(Γ)
    S₀ₕ = DiscreteSIO(S₀; h_mesh=0.1, h_quad=0.05)
    ϕ₀ₕ = S₀ₕ\f
    return true
end

function screen_test(Γ::SelfSimilarFractal)
    k = 1 + rand()*10
    f(x) = sin(norm(x)^2) # Laplace data

    amb_dim = Γ.spatial_dimension+1
    d = rand(amb_dim) # random 
    d ./ norm(d) # incident direction, normalised
    if amb_dim == 2
        θ = rand()
        x = [Γ.barycentre, 0] + Γ.diameter*(ones(amb_dim)+rand(amb_dim)) # random test point
    elseif amb_dim == 3
        θ = rand(2)
        x = [Vector{Float64}(Γ.barycentre); 0] + Γ.diameter*(ones(amb_dim)+rand(amb_dim)) # random test point
    end

    if amb_dim == 3
        g_ = (x::Vector{Float64}) -> exp(im*k*(d[1:(amb_dim-1)]'*x)) # screen pw data
    elseif amb_dim == 2
        g_ = (x::Float64) -> exp(im*k*d[1]*x) # screen pw data
    end

    Sₖ_ = SingleLayerOperatorHelmholtz(Γ, k; ambient_dimension=amb_dim)
    Sₖₕ_ = DiscreteSIO(Sₖ_; h_mesh=0.1, h_quad=0.05)
    ϕₖₕ_ = Sₖₕ_\g_
    𝙎ₖϕₕ = SingleLayerPotentialHelmholtz(ϕₖₕ_,k; ambient_dimension=amb_dim, h_quad = 0.5)
    𝙎ₖϕₕ(x)

    Sₖ_ = SingleLayerOperatorHelmholtz(Γ, k; ambient_dimension=amb_dim)
    Sₖₕ_ = DiscreteSIO(Sₖ_; h_mesh=0.1, h_quad=0.05)
    ϕₖₕ_ = Sₖₕ_\g_
    𝙁ₖϕₕ = far_field_pattern(ϕₖₕ_,k; h_quad = 0.05)
    𝙁ₖϕₕ(θ)

    S₀_ = SingleLayerOperatorLaplace(Γ; ambient_dimension=amb_dim)
    S₀ₕ_ = DiscreteSIO(S₀_;h_mesh=0.1, h_quad=0.05)
    ϕ₀ₕ_ = S₀ₕ_\f
    return true
end

function VIE_test(Γ::SelfSimilarFractal)
    k = 1 + rand()*10

    amb_dim = Γ.spatial_dimension
    x = Vector{Float64}(Γ.barycentre) + Γ.diameter*(ones(amb_dim)+rand(amb_dim)) # random test point
    d = rand(amb_dim) # random 
    d ./ norm(d) # incident direction, normalised

    nᵢ = 1+rand()
    uᴵ(x) = exp(im*k*(d'*x))

    VIO = Id(Γ) + k^2*(1-nᵢ)*VolumePotential(Γ,k)
    VIOₕ = DiscreteSIO(VIO, h_mesh =0.1, h_quad = 0.05)
    uᵀₕ = VIOₕ \ uᴵ
    Vuᵀₕ = IFSintegrals.volume_potential(uᵀₕ, k, h_quad=0.05)
    uᵀₕ_fn(x::Vector{Float64}) = uᴵ(x)-k^2*(1-nᵢ)*Vuᵀₕ(x)
    uᵀₕ_fn(x)
    return true
end