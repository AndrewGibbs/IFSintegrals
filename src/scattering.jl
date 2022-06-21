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
    Y, W = full_surface_quad(density,h_quad)
    # if ϕₕ.domain.topological_dimension == 2
        function R(x_::Vector{Float64})
            R_ = zeros(Float64,length(Y::Vector{Vector{Float64}}))
            for n=1:length((Y::Vector{Vector{Float64}}))
                R_[n] = norm([x_[1]-(Y::Vector{Vector{Float64}})[n][1],x_[2]-(Y::Vector{Vector{Float64}})[n][2],x_[3]])
            end
            return R_
        end
        Sϕ(x::Vector{Float64}) = (W::Vector{Complex{Float64}})'*(exp.(im*k*R(x))./(4π*(R(x))))
    # elseif ϕₕ.domain.topological_dimension == 1
    #     function R(x_::Vector{Float64})
    #         R_ = zeros(Float64,length(Y::Vector{Float64}))
    #         for n=1:length((Y::Vector{Float64}))
    #             R_[n] = norm([x_[1]-(Y::Vector{Float64})[n][1],x_[2]])
    #         end
    #         return R_
    #     end
    #     Sϕ(x::Vector{Float64}) = (W::Vector{Complex{Float64}})'*im/4*besselh(0,1,k*R(x))
    # end
    return Sϕ
end

# simplify for Laplace kernel case:
single_layer_potential(Γ::SelfSimilarFractal,density::Projection,points::Array{<:Real,2}) = single_layer_potential(Γ,0.0,density,points)

# far-field kernel 2D:
FF2D_kernel(θ::Float64, x::Float64,k::Float64) = exp(-im*k*(cos(θ)*x))

# far-field kernel 3D:
FF3D_kernel(θ::Float64, ψ::Float64, x::Vector{Float64}, k::Float64) = exp.(-im*k*(cos(θ)*x[1]+sin(ψ)*x[2]))

function far_field_pattern(ϕₕ::Projection, k::Real=0.0; h_quad::Float64=0.1/k)
    Y, W = full_surface_quad(ϕₕ,h_quad)
    if ϕₕ.domain.topological_dimension == 1
        Fϕ = (θ::Float64) -> transpose(W::Vector{Complex{Float64}})*FF2D_kernel.(θ,Y,k)
    elseif ϕₕ.domain.topological_dimension == 2
        Fϕ = (θ::Vector{Float64}) -> transpose(W::Vector{Complex{Float64}})*FF3D_kernel.(θ[1],θ[2],Y,k)
    else
        error("Cannot compute far-field for this number of spatial dimensions")
    end
    return Fϕ
end

