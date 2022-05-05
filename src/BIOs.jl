# import Base: -
# import Base: \ # to be overloaded with discrete operators
# using ProgressMeter

struct BIO{V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    domain::SelfSimilarFractal{V,M}
    kernel::Function #function
    Lipschitz_part_of_kernel::Function #function
    singularity_strength::Real
    singularity_scale::Complex{<:Real}
    self_adjoint::Bool
    coercive:: Bool
    wavenumber::Number
end

#constructor for zero wavenumber case
BIO(domain::SelfSimilarFractal{V,M},kernel::Function,Lipschitz_part_of_kernel::Function,singularity_strength::Real,
singularity_scale::Complex{<:Real},self_adjoint::Bool,coercive::Bool) where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}} =
BIO{V,M}(domain,kernel,Lipschitz_part_of_kernel,singularity_strength,singularity_scale,self_adjoint,coercive,0.0)

"""
    DiscreteBIO(BIO::BIO, h_BEM::Real, h_quad::Real)
    
is the constructor for a discretisation of a boundary layer boundary integral operator, 'BIO'.
h_BEM is the meshwidth parameter for the discretisation of the underlying fractal
h_quad denotes the discretisation parameter for the integrals in the stiffness matrix.
"""
struct DiscreteBIO{V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    BIO::BIO{V,M}
    h_BEM::Float64
    h_quad::Float64
    mesh::Vector{SubAttractor{V,M}}
    Lₕ::Vector{Vector{Int64}} # subindices list
    Galerkin_matrix::Matrix{<:Complex{<:Float64}} # not sure how to parametrise this as subtype of Array{<:Complex{<:Real},2}
end

#constructor:
function DiscreteBIO(K::BIO; h_BEM::Real=max(2π/(10.0*K.wavenumber),K.domain.diameter+eps()), h_quad::Real=h_BEM, Cosc = 2π, vary_quad=true, repeat_blocks=true)
    Γ = K.domain
    Lₕ = subdivide_indices(K.domain,h_BEM)
    N = length(Lₕ)
    mesh = [SubAttractor(Γ,Lₕ[n]) for n=1:N]
    M = length(Γ.IFS)
    #initialise Galerkin matrix:
    Galerkin_matrix = zeros(Complex{Float64},N,N)
    # create blank matrix of flags, describing if the matrix entry has been filled
    BEM_filled = zeros(Bool,N,N)
    m_count = 0

    if vary_quad
        h_quad_adjust = get_quad_scales(K,Lₕ)
    else
        h_quad_adjust = ones(Float64, length(Lₕ),length(Lₕ))
    end

    if Γ.homogeneous & repeat_blocks
        ℓ = length(Lₕ[1])
        diag_block_sizes = M.^(0:(ℓ-1))
    else
        diag_block_sizes = []
    end

    @showprogress 1 "Constructing BEM system " for m_count=1:N#m in Lₕ
        m = Lₕ[m_count]
        Γₘ = mesh[m_count]

        # if matrix is symmetric, will only need to compute ≈ half entries
        K.self_adjoint ? n_count_start = m_count : n_count_start = 1

        for n_count = n_count_start:length(Lₕ)#n in Lₕ
            if !BEM_filled[m_count,n_count] # check matrix entry hasn't been filled already
                n = Lₕ[n_count]
                Γₙ = mesh[n_count]
                x,y,w = barycentre_rule(Γₘ,Γₙ,h_quad*h_quad_adjust[m_count,n_count])
                if n==m
                    Galerkin_matrix[m_count,n_count] = singular_elliptic_double_integral(K,h_quad,n;Cosc=Cosc)
                else
                    Galerkin_matrix[m_count,n_count] = w'*K.kernel.(x,y)
                end

                # if matrix is symmetric, expoit this to save time
                if K.self_adjoint && n!=m
                    Galerkin_matrix[n_count,m_count] = Galerkin_matrix[m_count,n_count]
                end
            end
        end

        # now repeat entries along block diagonal, if at the end of a diagonal block, and homogeneous
        if in(m_count,diag_block_sizes) && N>1
            block_power = indexin(m_count,diag_block_sizes)[1]
            block_size = M^(block_power-1)
            for j=2:M
                block_range = ((j-1)*block_size+1):(j*block_size)
                Galerkin_matrix[block_range,block_range] = Galerkin_matrix[1:block_size,1:block_size]
                BEM_filled[block_range,block_range] .= 1
            end
        end
    end
    DiscreteBIO(K, h_BEM, h_quad, mesh, Lₕ, Galerkin_matrix)
end

"""
    SingleLayer(Γ::SelfSimilarFractal, k::Real=0.0)

represents the single layer boundary integral operator, Sϕ(x) = ∫_Γ Φ(x,y) ϕ(x) dHᵈ(y),
where Φ is the fundamental solution for the underlying PDE.
"""
function SingleLayer(Γ::SelfSimilarFractal{V,M}, k::Number=0.0) where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    if Γ.topological_dimension == 1
        if k==0.0 #2D Laplace case
            K = BIO{V,M}(Γ, #fractal domain
            (x,y)->Φₜ(0.0,x,y), # log kernel
            (x,y)->zero_kernel(x,y), # kernel minus singularity
            0.0, # strength of singularity, corresponding to log singularity
            -1/(2π), # scaling of singularity
            true, #self-adjoint
            true, #coercive
            k #wavenumber
            )
        else #2D Helmholtz case        
            K = BIO{V,M}(Γ, #fractal domain
            (x,y)->HelhmoltzGreen2D(k,x,y), # Hankel function
            (x,y)->HelhmoltzGreen2D_Lipschitz_part(k,x,y), # kernel minus singularity
            0.0, # strength of singularity, corresponding to log singularity
            -1/(2π), # scaling of singularity
            true, #self-adjoint
            true, #coercive
            k #wavenumber
            )
        end
    elseif Γ.topological_dimension == 2
        if k==0.0 #3D Laplace case
            K = BIO{V,M}(Γ, #fractal domain
            (x,y)-> Φₜ(1.0,x,y), # Green's function
            (x,y)-> zero_kernel(x,y), # kernel minus singularity
            1.0, # strength of singularity, corresponding to 1/|x-y|
            1/(4π), # scaling of singularity
            true, #self-adjoint
            true, #coercive
            k #wavenumber
            )
        else #3D Helmholtz case        
            K = BIO{V,M}(Γ, #fractal domain
            (x,y)->HelhmoltzGreen3D(k,x,y), # Green's function
            (x,y)->HelhmoltzGreen3D_Lipschitz_part(k,x,y), # kernel minus singularity
            1.0, # strength of singularity, corresponding to 1/|x-y|
            1/(4π), # scaling of singularity
            true, #self-adjoint
            true, #coercive
            k #wavenumber
            )
        end
    else
        error("Haven't coded single layer BIO for this many dimensions")
    end
end

function singular_elliptic_double_integral(K::BIO, h_quad::Real,index::Array{Int64}=[0]; Cosc = 2π)
    if real(K.wavenumber)*SubAttractor(K.domain,index).diameter > Cosc
        return singular_elliptic_double_integral_full(K::BIO,h_quad::Real,index; Cosc = Cosc)
    else
        return singular_elliptic_double_integral_basic(K::BIO,h_quad::Real,index)
    end
end

function singular_elliptic_double_integral(Γ::Union{Attractor,SubAttractor},k::Number,h_quad::Real; Cosc = 2π)
    K = SingleLayer(Γ, k)
    return singular_elliptic_double_integral(K, h_quad; Cosc = Cosc)
end

function singular_elliptic_double_integral_full(K::BIO,h_quad::Real,index::Array{Int64}=[0]; Cosc = 2π)
    # following the procedure in the paper, to avoid errors in far-field.
    # notation agrees with quadrature paper here, rather than BEM paper.
    Γ = SubAttractor(K.domain,index)
    h_star = Cosc/real(K.wavenumber)
    h = min(h_star,h_quad)
    L_h_star = subdivide_indices(Γ,h_star)
    Npts = 0

    I = zero(Complex{Float64})
    for n in L_h_star
        Γₙ = SubAttractor(Γ,n)
        for m in L_h_star
            Γₘ = SubAttractor(Γ,m)
            x,y,w = barycentre_rule(Γₙ,Γₘ,h)
            Npts += length(w) # these points will get used at 
            if m == n # diagonal quadrature element
                I += K.singularity_scale*eval_green_double_integral(Γₘ,K.singularity_strength,h)
                I += w'*K.Lipschitz_part_of_kernel.(x,y)
            else
                I += w'*K.kernel.(x,y)
            end
        end
    end
    return I
end

function singular_elliptic_double_integral_basic(K::BIO,h::Real,index::Array{Int64}=[0])
    Γ = SubAttractor(K.domain,index)
    x,y,w = barycentre_rule(Γ,Γ,h)
    I = K.singularity_scale*eval_green_double_integral(Γ,K.singularity_strength,h) + w'*K.Lipschitz_part_of_kernel.(x,y)
    return I
end

# The next couple of functions are designed to use less quadrature points in the BEM elements
# where this won't affect the accuracy.

F_nomeasure(r::Real, k::Number, n::Int64) = (1+(abs(k)*r)^(n/2+1))/r^(n+1)

function get_quad_scales(K::BIO,Lₕ::Vector{Vector{Int64}})
    Γ = K.domain
    # compute upper and lower bounds for the F in my notes, which is stated above.
    F_upper = ones(Float64,length(Lₕ),length(Lₕ))
    F_lower = ones(Float64,length(Lₕ),length(Lₕ))
    for m_count = 1:length(Lₕ)
        Γₘ = SubAttractor(Γ,Lₕ[m_count])
        for n_count = m_count:length(Lₕ)
            Γₙ = SubAttractor(Γ,Lₕ[n_count])
            if n_count!=m_count
                dist_upper = norm(Γₙ.barycentre-Γₘ.barycentre)
                dist_lower = max(dist_upper - Γₙ.diameter - Γₘ.diameter,0)
                measure_weight = Γₙ.measure*Γₘ.measure
                # noting that F_nomeasure is monotonic decreasing in r, can bound as follows:
                if dist_lower>0
                    F_upper[m_count,n_count] = measure_weight*F_nomeasure(dist_lower, K.wavenumber, K.domain.topological_dimension)
                else
                    F_upper[m_count,n_count]= Inf
                end
                 F_lower[m_count,n_count] = measure_weight*F_nomeasure(dist_upper, K.wavenumber, K.domain.topological_dimension)
            else
                F_upper[m_count,n_count] = Inf
            end
            
            # now by symmetry
            F_upper[n_count,m_count] = F_upper[m_count,n_count]
            F_lower[n_count,m_count] = F_lower[m_count,n_count]
        end
    end
    # now can get a lower bound estimate on the quantity from my notes:
    quad_scales = max.(sqrt.(maximum(F_lower)./F_upper),1)
    return quad_scales
end
