import Base:\ # to be overloaded with discrete operators
using ProgressMeter

struct BIO
    domain::Fractal
    kernel #function
    Lipschitz_part_of_kernel #function
    singularity_strength::Real
    singularity_scale::Complex{<:Real}
    wavenumber::Real
end

"""
    DiscreteBIO(BIO::BIO, h_BEM::Real, h_quad::Real)
    
is the constructor for a discretisation of a boundary layer boundary integral operator, 'BIO'.
h_BEM is the meshwidth parameter for the discretisation of the underlying fractal
h_quad denotes the discretisation parameter for the integrals in the stiffness matrix.
"""
struct DiscreteBIO
    BIO::BIO
    h_BEM::Real
    h_quad::Real
    Lₕ # subindices list
    Galerkin_matrix::Array{<:Complex{<:Real},2} # not sure how to parametrise this as subtype of Array{<:Complex{<:Real},2}
end

#constructor:
function DiscreteBIO(K::BIO, h_BEM::Real, h_quad::Real)
    Γ = K.domain
    Lₕ = subdivide_indices(K.domain,h_BEM)
    N = length(Lₕ)
    #initialise Galerkin matrix:
    Galerkin_matrix = zeros(Complex{Float64},N,N)
    m_count = 0
    @showprogress 1 "Constructing BEM system " for m in Lₕ
        m_count+=1
        Γₘ = SubAttractor(Γ,m)
        n_count = 0
        for n in Lₕ
            n_count += 1
            Γₙ = SubAttractor(Γ,n)
            x,y,w = barycentre_rule(Γₘ,Γₙ,h_quad)
            if n==m
                Galerkin_matrix[m_count,n_count] = singular_elliptic_double_integral(K,h_quad,n)
            else
                Galerkin_matrix[m_count,n_count] = sum(w.*K.kernel(x,y))
            end
        end
    end
    DiscreteBIO(K, h_BEM, h_quad, Lₕ, Galerkin_matrix)
end

function \(K::DiscreteBIO, f::Function)
    Γ = K.BIO.domain
    F = zeros(Complex{Float64},length(K.Lₕ))
    n_count = 0
    for n in K.Lₕ
        n_count += 1
        Γₙ = SubAttractor(Γ, n)
        x,w = barycentre_rule(Γₙ, K.h_quad)
        X = slicematrix(x)
        F[n_count] = sum(w.*f.(X))
    end
    coeffs = K.Galerkin_matrix \ F
    return Projection(K.BIO.domain,K.Lₕ,coeffs)
end

struct Projection # onto the coefficient space of piecewise constants on the fractal
    domain::Fractal
    Lₕ # subindices list
    coeffs::Array{<:Complex{<:Real},1}
end

"""
    SingleLayer(Γ::Fractal, k::Real=0.0)

represents the single layer boundary integral operator, Sϕ(x) = ∫_Γ Φ(x,y) ϕ(x) dHᵈ(y),
where Φ is the fundamental solution for the underlying PDE.
"""
function SingleLayer(Γ::Fractal, k::Real=0.0)
    if Γ.topological_dimension == 1
        if k==0.0 #2D Laplace case
            K = BIO(Γ, #fractal domain
            (x,y)->Φₜ(0.0,x,y), # log kernel
            (x,y)->zero_kernel(x,y), # kernel minus singularity
            0.0, # strength of singularity, corresponding to log singularity
            -1/(2π), # scaling of singularity
            k #wavenumber
            )
        else #2D Helmholtz case        
            K = BIO(Γ, #fractal domain
            (x,y)->HelhmoltzGreen2D(k,x,y), # Hankel function
            (x,y)->HelhmoltzGreen2D_Lipschitz_part(k,x,y), # kernel minus singularity
            0.0, # strength of singularity, corresponding to log singularity
            -1/(2π), # scaling of singularity
            k #wavenumber
            )
        end
    elseif Γ.topological_dimension == 2
        if k==0.0 #3D Laplace case
            K = BIO(Γ, #fractal domain
            (x,y)-> Φₜ(1.0,x,y), # Green's function
            (x,y)-> zero_kernel(x,y), # kernel minus singularity
            1.0, # strength of singularity, corresponding to 1/|x-y|
            1/(4π), # scaling of singularity
            k #wavenumber
            )
        else #3D Helmholtz case        
            K = BIO(Γ, #fractal domain
            (x,y)->HelhmoltzGreen3D(k,x,y), # Green's function
            (x,y)->HelhmoltzGreen3D_Lipschitz_part(k,x,y), # kernel minus singularity
            1.0, # strength of singularity, corresponding to 1/|x-y|
            1/(4π), # scaling of singularity
            k #wavenumber
            )
        end
    else
        error("Haven't coded single layer BIO for this many dimensions")
    end
end

function singular_elliptic_double_integral(K::BIO,h_quad::Real,index::Array{Int64}=[0])
    # following the procedure in the paper, to avoid errors in far-field.
    # notation agrees with quadrature paper here, rather than BEM paper.
    # if index == [0]
    #     Γ = K.domain
    # else
    #     Γ = SubAttractor(K.domain,index)
    # end
    # now partition into subintervals, based on this frequency-dependent parameter:
    Γ = SubAttractor(K.domain,index)
    h_star = min(1/K.wavenumber,Γ.diameter)
    h = min(h_star,h_quad)
    L_h_star = subdivide_indices(Γ,h_star)

    I = zero(Complex{Float64})
    for n in L_h_star
        Γₙ = SubAttractor(Γ,n)
        for m in L_h_star
            Γₘ = SubAttractor(Γ,m)
            x,y,w = barycentre_rule(Γₙ,Γₘ,h)
            if m == n # diagonal quadrature element
                I += K.singularity_scale*eval_green_double_integral(Γₘ,K.singularity_strength,h)
                I += sum(w.*K.Lipschitz_part_of_kernel(x,y))
            else
                I += sum(w.*K.kernel(x,y))
            end
        end
    end
    return I
end

function singular_elliptic_double_integral_old(K::BIO,h::Real,index::Array{Int64}=[0])
    # following the procedure in the paper, to avoid errors in far-field.
    # notation agrees with quadrature paper here, rather than BEM paper.
    # if index == [0]
    #     Γ = K.domain
    # else
    #     Γ = SubAttractor(K.domain,index)
    # end
    # now partition into subintervals, based on this frequency-dependent parameter:
    Γ = SubAttractor(K.domain,index)
    x,y,w = barycentre_rule(Γ,Γ,h)
    I = K.singularity_scale*eval_green_double_integral(Γ,K.singularity_strength,h) + sum(w.*K.Lipschitz_part_of_kernel(x,y))
    return I
end

# function single_layer_potential(Γ::Fractal,k::Real,density::Projection,points::Array{<:Array{<:Real}})

# end

function single_layer_potential(Γ::Fractal,k::Real,density::Projection,points::Array{<:Real,2})
    h = 0.1/k
    num_points = size(points)[1]
    vals = zeros(Complex{Float64},num_points)
    # kernel is the same as the single layer BIO, so use that:
    S = SingleLayer(Γ, k)
    Φ = S.kernel
    for m_count =1:length(density.Lₕ)
        m = density.Lₕ[m_count]
        Γₘ = SubAttractor(Γ,m)
        y,w = barycentre_rule(Γₘ,h)
        # now convert y to n+1 dimensions
        Y = [y zeros(length(w))]
        for n = 1:num_points
            vals[n] += w'*Φ(Y,points[n,:])*density.coeffs[m_count]
        end
        m_count += 1
    end
    return vals
end
# simplify for Laplace kernel case:
single_layer_potential(Γ::Fractal,density::Projection,points::Array{<:Real,2}) = single_layer_potential(Γ,0.0,density,points)