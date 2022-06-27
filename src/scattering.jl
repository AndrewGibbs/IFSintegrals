function full_surface_quad_N1(density::Projection,h::Float64)
    W = ComplexF64[]
    Y = Float64[]

    for m_count = 1:length(density.mesh)
        y,w = barycentre_rule(density.mesh[m_count],h)
        Y = vcat(Y,y)
        # z = Vector{ComplexF64}(w*density.coeffs[m_count])
        W = vcat(W,w*density.coeffs[m_count])
    end
    return Y,W
end

function full_surface_quad(density::Projection,h::Float64)
    W = Complex{Float64}[]
    Y = Vector{Float64}[]

    for m_count = 1:length(density.mesh)
        y,w = barycentre_rule(density.mesh[m_count],h)
        Y = vcat(Y,y)
        W = vcat(W,w*density.coeffs[m_count])
    end
    return Y,W
end

function single_layer_potential(density::Projection, k::Real; h_quad::Float64=0.1/k)
    # get points on surface of Γ
    if density.domain.topological_dimension == 1
        Y, W = full_surface_quad_N1(density,h_quad)
    else
        Y, W = full_surface_quad(density,h_quad)
    end
    # construct function for distance from screen
    function R(x_::Vector{Float64})
        R_ = (x_[end]^2)*ones(Float64,length(W))
        for n=1:length(W)
            for j = 1:(length(x_)-1)
                R_[n] += (x_[j]-Y[n][j]).^2
            end
        end
        return sqrt.(R_)
    end
    # create map, based on appropriate kernel
    if density.domain.topological_dimension == 1
        Sϕ = (x::Vector{Float64}) -> (W::Vector{Complex{Float64}})'*im/4*besselh.(0,1,k*R(x))
    elseif density.domain.topological_dimension == 2
        Sϕ = (x::Vector{Float64}) -> (W::Vector{Complex{Float64}})'*(exp.(im*k*R(x))./(4π*(R(x))))
    else
        error("Cannot compute single layer potential for this number of spatial dimensions")
    end
    return Sϕ
end

# far-field kernel 2D:
FF2D_kernel(θ::Float64, x::Float64,k::Float64) = exp(-im*k*(cos(θ)*x))

# far-field kernel 3D:
FF3D_kernel(θ::Float64, ψ::Float64, x::Vector{Float64}, k::Float64) = exp.(-im*k*(cos(θ)*x[1]+sin(ψ)*x[2]))

function far_field_pattern(ϕₕ::Projection, k::Real=0.0; h_quad::Float64=0.1/k)
    # get points on surface of Γ
    if ϕₕ.domain.topological_dimension == 1
        Y, W = full_surface_quad_N1(density,h_quad)
    else
        Y, W = full_surface_quad(density,h_quad)
    end
    # create map, based on appropriate kernel
    if ϕₕ.domain.topological_dimension == 1
        Fϕ = (θ::Float64) -> transpose(W::Vector{Complex{Float64}})*FF2D_kernel.(θ,Y,k)
    elseif ϕₕ.domain.topological_dimension == 2
        Fϕ = (θ::Vector{Float64}) -> transpose(W::Vector{Complex{Float64}})*FF3D_kernel.(θ[1],θ[2],Y,k)
    else
        error("Cannot compute far-field for this number of spatial dimensions")
    end
    return Fϕ
end

