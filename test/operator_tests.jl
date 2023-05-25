using IFSintegrals, LinearAlgebra, Test

macro no_error(ex)
    quote
        try
            $(esc(ex))
            true
        catch
            false
        end
    end
end

k = 1 + rand()*10

f(x) = sin(norm(x)^2) # Laplace data

screens = [CantorSet(), CantorDust(), Sierpinski(), KochFlake()]
surfaces = [CantorDust()]
volumes = [KochFlake(), Dragon()]

@testset "operator and potential tests" begin
    for Γ = surfaces

        amb_dim = Γ.spatial_dimension
        x = Γ.barycentre + Γ.diameter*(ones(amb_dim)+rand(amb_dim)) # random test point
        d = rand(amb_dim) # random 
        d ./ norm(d) # incident direction, normalised
        if amb_dim == 2
            θ = rand()
        elseif amb_dim == 2
            θ = rand(2)
        end

        g(x) = exp(im*k*(d'*x)) # pw data

        @test@no_error(Sₖ = SingleLayerOperatorHelmholtz(Γ, k))
        @test@no_error(Sₖₕ = DiscreteSIO(Sₖ; h_mesh=0.1, h_quad=0.05))
        @test@no_error(ϕₖₕ = Sₖₕ\g)
        @test@no_error(𝙎ₖϕₕ = SingleLayerPotentialHelmholtz(ϕₕ,k; h_quad = 0.05))
        @test@no_error(𝙎ₖϕₕ(x))

        @test@no_error(𝙁ₖϕₕ = far_field_pattern(ϕₕ,k; h_quad = 0.05))
        @test@no_error(𝙁ₖϕₕ(θ))

        @test@no_error(S₀ = SingleLayerOperatorLaplace(Γ, k))
        @test@no_error(S₀ₕ = DiscreteSIO(Sₖ;h_mesh=0.1, h_quad=0.05))
        @test@no_error(ϕ₀ₕ = Sₖₕ\f)

    end
# end

# @testset "screen BEMs" begin
    for Γ = screens

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

        g_(x) = exp(im*k*(d[1:(amb_dim-1)]'*x)) # screen pw data

        @test@no_error(Sₖ_ = SingleLayerOperatorHelmholtz(Γ, k; ambient_dimension=amb_dim))
        @test@no_error(Sₖₕ_ = DiscreteSIO(Sₖ; h_mesh=0.1, h_quad=0.05))
        @test@no_error(ϕₖₕ_ = Sₖₕ\g_)
        @test@no_error(𝙎ₖϕₕ = SingleLayerPotentialHelmholtz(ϕₕ,k; ambient_dimension=amb_dim, h_quad = 0.05))
        @test@no_error(𝙎ₖϕₕ(x))

        @test@no_error(𝙁ₖϕₕ = far_field_pattern(ϕₕ,k; h_quad = 0.05))
        @test@no_error(𝙁ₖϕₕ(θ))

        @test@no_error(S₀_ = SingleLayerOperatorLaplace(Γ, k; ambient_dimension=amb_dim))
        @test@no_error(S₀ₕ_ = DiscreteSIO(Sₖ;h_mesh=0.1, h_quad=0.05))
        @test@no_error(ϕ₀ₕ_ = Sₖₕ\f)

    end
# end

# @testset "VIEs" begin
    for Ω = volumes

        amb_dim = Γ.spatial_dimension
        x = Γ.barycentre + Γ.diameter*(ones(amb_dim)+rand(amb_dim)) # random test point
        d = rand(amb_dim) # random 
        d ./ norm(d) # incident direction, normalised

        nᵢ = 1+rand()
        uᴵ(x) = exp(im*k*(d'*x))

        @test@no_error(VIO = Id(Γ) + k^2*(1-nᵢ)*VolumePotential(Γ,k))
        @test@no_error(VIOₕ = DiscreteSIO(VIO, h_mesh = h_mesh, h_quad = h_quad))
        @test@no_error(uᵀₕ = VIOₕ \ uᴵ)
        @test@no_error(Vuᵀₕ = IFSintegrals.volume_potential(uᵀₕ, k, h_quad=h_quad))
        @test@no_error(uᵀₕ_fn(x::Vector{Float64}) = uᴵ(x)-k^2*(1-nᵢ)*Vuᵀₕ(x))
        @test@no_error(uᵀₕ_fn(x))
    end
end