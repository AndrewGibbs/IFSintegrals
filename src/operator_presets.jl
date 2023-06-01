HelhmoltzGreen2D(k::Number,x::T, y::T) where T<:Union{Real,AbstractVector} = im/4*besselh(0,1,k*norm(x-y))
HelhmoltzGreen2D(k::Number,r::Real) = im/4*besselh(0,1,k*r)

function HelhmoltzGreen2D_Lipschitz_part(k::Number, x::T, y::T) where T<:Union{Real,AbstractVector}
    if x == y
        return im/4 -1/(2π)*(0.577215664901532 + log(k/2))
    else
        return HelhmoltzGreen2D(k,x,y) + 1/(2π)*log(norm(x-y))
    end
end

HelhmoltzGreen3D(k::Number,x::T,y::T) where T<:Union{Real,AbstractVector} = exp(im*k*norm(x-y))/(4π*norm(x-y))
HelhmoltzGreen3D(k::Number,r::Real) = exp(im*k*r)/(4π*r)

function HelhmoltzGreen3D_Lipschitz_part(k::Number, x::T, y::T) where T<:Union{Real,AbstractVector}
    if x == y
        return im*k/(4π)
    else
        return expm1(im*k*norm(x-y)) /(4π*norm(x-y))#HelhmoltzGreen3D(k,x,y) - 1.0 /(4π*norm(x-y))
    end
end

# get_spatial_dimension(Γ::SelfSimilarFractal)
"""
    SingleLayerOperatorLaplace(Γ::SelfSimilarFractal, wavenumber::Real=0.0)

represents the single layer boundary integral operator for Laplace, Sϕ(x) = ∫_Γ Φ(x,y) ϕ(x) dHᵈ(y),
where Φ is the fundamental solution for the underlying PDE.
"""
function SingleLayerOperatorLaplace(Γ::Ω; ambient_dimension::Int64 = Γ.spatial_dimension
        ) where Ω <: SelfSimilarFractal# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    if ambient_dimension == 2
        K = SIO{Ω}(Γ, #fractal domain
        (x,y)->Φₜ(0.0,x,y), # log kernel
        (x,y)->zero_kernel(x,y), # kernel minus singularity
        0.0, # strength of singularity, corresponding to log singularity
        -1/(2π), # scaling of singularity
        true, #self-adjoint
        0.0 #wavenumber
        )
    elseif ambient_dimension == 3
        K = SIO{typeof(Γ)}(Γ, #fractal domain
        (x,y)-> Φₜ(0.0,x,y), # Green's function
        (x,y)-> zero_kernel(x,y), # kernel minus singularity
        1.0, # strength of singularity, corresponding to 1/|x-y|
        1/(4π), # scaling of singularity
        true, #self-adjoint
        0.0 #wavenumber
        )
    else
        error("Haven't coded single layer SIO for this many dimensions")
    end
end

"""
    SingleLayerOperatorHelmholtz(Γ::SelfSimilarFractal, wavenumber::Real=0.0)

represents the single layer boundary integral operator for Helmholtz, Sϕ(x) = ∫_Γ Φ(x,y) ϕ(x) dHᵈ(y),
where Φ is the fundamental solution for the underlying PDE.
"""
function SingleLayerOperatorHelmholtz(Γ::Ω, k::Number; ambient_dimension::Int64 = Γ.spatial_dimension
    ) where Ω <: SelfSimilarFractal
    if ambient_dimension == 2     
        K = SIO{Ω}(Γ, #fractal domain
        (x,y)->HelhmoltzGreen2D(k,x,y), # Hankel function
        (x,y)->HelhmoltzGreen2D_Lipschitz_part(k,x,y), # kernel minus singularity
        0.0, # strength of singularity, corresponding to log singularity
        -1/(2π), # scaling of singularity
        true, #self-adjoint
        k #wavenumber
        )
    elseif ambient_dimension == 3
        #3D Helmholtz case        
            K = SIO(Γ, #fractal domain
            (x,y)->HelhmoltzGreen3D(k,x,y), # Green's function
            (x,y)->HelhmoltzGreen3D_Lipschitz_part(k,x,y), # kernel minus singularity
            1.0, # strength of singularity, corresponding to 1/|x-y|
            1/(4π), # scaling of singularity
            true, #self-adjoint
            k #wavenumber
            )
    else
        error("Haven't coded single layer SIO for this many dimensions")
    end
end


"""
    single_layer_potential(density::Projection, wavenumber::Real; h_quad::Float64=0.1/k)
Returns a function which is the single layer potential in Rⁿ⁺¹, from a density on the screen Γ ∈ Rⁿ.
"""
function get_layer_potential(density::Projection, Φᵣ::Function, h_quad)# where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    # get points on surface of Γ
    Y, W = full_surface_quad(density,h_quad)
    vec_length = length(W)
    spatial_dims = density.domain.spatial_dimension

    function Sϕ(x::Vector{Float64})
        R = zeros(Float64, vec_length)
        # first get contributions from additional spatial dimensions, if any
        for j=(spatial_dims+1):length(x)
            R .+=x[j]^2
        end
        for n=1:vec_length
            for j = 1:spatial_dims
                R[n] += (x[j]-Y[n][j]).^2
            end
        end
        R = Vector{Float64}(sqrt.(R))
        return transpose(W::Vector{Complex{Float64}})*(Φᵣ.(R))
    end
    return Sϕ
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