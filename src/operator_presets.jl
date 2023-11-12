# get_spatial_dimension(Γ::FractalMeasure)
"""
    SingleLayerOperatorLaplace(Γ::FractalMeasure, wavenumber::Real=0.0)

represents the single layer boundary integral operator for Laplace, Sϕ(x) = ∫_Γ Φ(x,y) ϕ(x) dHᵈ(y),
where Φ is the fundamental solution for the underlying PDE.
"""
function SingleLayerOperatorLaplace(Γ::Ω; ambient_dimension::Int64 = Γ.spatial_dimension
        ) where Ω <: FractalMeasure# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    if ambient_dimension == 2
        K = SIO{Ω}(Γ, #fractal domain
        (x,y)->-1/(2π)*Φₜ(0.0,x,y), # log kernel
        (x,y)->zero_kernel(x,y), # kernel minus singularity
        0.0, # strength of singularity, corresponding to log singularity
        -1/(2π), # scaling of singularity
        true, #self-adjoint
        0.0 #wavenumber
        )
    elseif ambient_dimension == 3
        K = SIO{Ω}(Γ, #fractal domain
        (x,y)-> 1/(4π)*Φₜ(1.0,x,y), # Green's function
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
    SingleLayerOperatorHelmholtz(Γ::FractalMeasure, wavenumber::Real=0.0)

represents the single layer boundary integral operator for Helmholtz, Sϕ(x) = ∫_Γ Φ(x,y) ϕ(x) dHᵈ(y),
where Φ is the fundamental solution for the underlying PDE.
"""
function SingleLayerOperatorHelmholtz(Γ::Ω, k::Number; ambient_dimension::Int64 = Γ.spatial_dimension
    ) where Ω <: FractalMeasure
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
            K = SIO{Ω}(Γ, #fractal domain
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

