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


"""
    single_layer_potential(density::Projection, wavenumber::Real; h_quad::Float64=0.1/k)
Returns a function which is the single layer potential in Rⁿ⁺¹, from a density on the screen Γ ∈ Rⁿ.
"""
# function single_layer_potential(density::Projection{V,M}, k::Real; h_quad::Float64=0.1/k) where {V<:AbstractVector, M<:Union{Real,AbstractMatrix}}
#     # get points on surface of Γ
#     Y, W = full_surface_quad(density,h_quad)
#     vec_length = length(W)
#     function Sϕ(x::Vector{Float64})
#         R = (x[3]^2)*ones(Float64,vec_length)
#         for n=1:vec_length
#             for j = 1:2
#                 R[n] += (x[j]-Y[n][j]).^2
#             end
#         end
#         R = Vector{Float64}(sqrt.(R))
#         return transpose(W::Vector{Complex{Float64}})*(exp.(im*k*R)./(4π*R))
#     end
#     return Sϕ
# end

# function single_layer_potential(density::Projection{V,M}, k::Real; h_quad::Float64=0.1/k) where {V<:Real, M<:Real}
#     density.domain.spatial_dimension > 2 ? error("Cannot compute single layer potential for this number of spatial dimensions") : Nothing
#     # get points on surface of Γ
#     Y, W = full_surface_quad_N1(density,h_quad)
#     vec_length = length(W)
#     function Sϕ(x::Vector{Float64})
#         R = (x[2]^2)*ones(Float64,vec_length)
#         for n=1:vec_length
#             R[n] += (x[1]-Y[n]).^2
#         end
#         R = Vector{Float64}(sqrt.(R))
#         return transpose(W::Vector{Complex{Float64}})*im/4*besselh.(0,1,k*R)
#     end
#     # error("Cannot compute single layer potential for this number of spatial dimensions")
#     return Sϕ
# end

# far-field kernel 2D:
FF2D_kernel_screen(θ::Float64, x::Float64,k::Number) = 
    -sqrt(1im/(8*π*k))*exp(-im*k*(cos(θ)*x))
FF2D_kernel(θ::Float64, x::Vector{Float64},k::Number) = 
    -sqrt(1im/(8*π*k))*exp(-im*k*(cos(θ)*x[1]+sin(θ)*x[2]))

# far-field kernel 3D:
FF3D_kernel_screen(θ::Float64, ψ::Float64, x::Vector{Float64}, k::Number) = 
    -(1/(4π))*exp(-im*k*(sin(θ)*cos(ψ)*x[1] + sin(θ)*sin(ψ)*x[2]))
FF3D_kernel(θ::Float64, ψ::Float64, x::Vector{Float64}, k::Number) = 
    -(1/(4π))*exp(-im*k*(sin(θ)*cos(ψ)*x[1] + sin(θ)*sin(ψ)*x[2] + cos(θ)*x[3]))

function VolumePotential(Γ::FractalMeasure, k::Number)# where Ω <: FractalMeasure
    @warn("VolumePotential will soon be depracted. Use SingleLayerOperatorHelmholtz instead.")
    return SingleLayerOperatorHelmholtz(Γ, k)
end

function SingleLayerPotentialHelmholtz(density::Projection,
    k::T; ambient_dimension::Int64 = density.domain.spatial_dimension,
    h_quad::Float64 = 0.1/abs(k)
    ) where T<:Number# where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}, T<:Number}

    # choose appropriate Green's kernel
    if ambient_dimension == 2
        K = (r::Float64)-> HelhmoltzGreen2D(k::T,r)
    elseif ambient_dimension == 3
        K = (r::Float64)-> HelhmoltzGreen3D(k::T,r)
    else
        error("Cannot compute single layer potential for this number of spatial dimensions")
    end
    # return the potential function
    return get_layer_potential(density, K, h_quad)
end

function far_field_pattern(density::Projection,
    k::T; ambient_dimension::Int64 = density.domain.spatial_dimension,
    h_quad::Float64 = 0.1/abs(k)
    ) where T<:Number

    if density.domain.spatial_dimension == 1
        Y, W = full_surface_quad_N1(density,h_quad)
    else
        Y, W = full_surface_quad(density,h_quad)
    end
    
    # choose appropriate Green's kernel
    if ambient_dimension == 2
        if density.domain.spatial_dimension == 1
            Fϕ = (θ::Float64) -> transpose(W::Vector{Complex{Float64}})*FF2D_kernel_screen.(θ,Y,k)
        else
            Fϕ = (θ::Float64) -> transpose(W::Vector{Complex{Float64}})*FF2D_kernel.(θ,Y,k)
        end
    elseif ambient_dimension == 3
        if density.domain.spatial_dimension == 2
            Fϕ = (θ::Vector{Float64}) -> transpose(W::Vector{Complex{Float64}})*FF3D_kernel_screen.(θ[1],θ[2],Y,k)
        else
            Fϕ = (θ::Vector{Float64}) -> transpose(W::Vector{Complex{Float64}})*FF3D_kernel.(θ[1],θ[2],Y,k)
        end
    else
        error("Cannot compute single layer potential for this number of spatial dimensions")
    end

    return Fϕ
end
"""
Returns a function which is the far-field pattern in Rⁿ⁺¹, from a density on the screen Γ ∈ Rⁿ.
"""
# function far_field_pattern_3D(density::Projection{V,M}, k::Real=0.0; h_quad::Float64=0.1/k) where {V<:AbstractVector, M<:Union{Real,AbstractMatrix}}
#     density.domain.spatial_dimension > 2 ? error("Cannot compute far-field pattern for this number of spatial dimensions") : Nothing# get points on surface of Γ
#     Y, W = full_surface_quad(density,h_quad)
#     # create map, based on appropriate kernel
#     Fϕ = (θ::Vector{Float64}) -> transpose(W::Vector{Complex{Float64}})*FF3D_kernel.(θ[1],θ[2],Y,k)
#     return Fϕ
# end

# function far_field_pattern_2D(density::Projection{V,M}, k::Real=0.0; h_quad::Float64=0.1/k) where {V<:Real, M<:Real}
#     # get points on surface of Γ
#     Y, W = full_surface_quad_N1(density,h_quad)
#     # create map, based on appropriate kernel
#     Fϕ = (θ::Float64) -> transpose(W::Vector{Complex{Float64}})*FF2D_kernel.(θ,Y,k)
#     return Fϕ
# end