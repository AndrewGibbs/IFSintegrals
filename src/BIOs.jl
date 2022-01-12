import Base: -
import Base: \ # to be overloaded with discrete operators
using ProgressMeter

struct BIO
    domain::Fractal
    kernel::Function #function
    Lipschitz_part_of_kernel::Function #function
    singularity_strength::Real
    singularity_scale::Complex{<:Real}
    self_adjoint::Bool
    wavenumber::Number
end

#constructor for zero wavenumber case
BIO(domain::Fractal,kernel::Function,Lipschitz_part_of_kernel::Function,singularity_strength::Real,
singularity_scale::Complex{<:Real},self_adjoint::Bool)=BIO(domain,kernel,Lipschitz_part_of_kernel,
singularity_strength,singularity_scale,self_adjoint,0.0)

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
    Lₕ::Array{Array{Int64,1},1} # subindices list
    Galerkin_matrix::Array{<:Complex{<:Real},2} # not sure how to parametrise this as subtype of Array{<:Complex{<:Real},2}
end

#constructor:
function DiscreteBIO(K::BIO, h_BEM::Real, h_quad::Real; Cosc = 2π, vary_quad=true, repeat_blocks=true)
    Γ = K.domain
    Lₕ = subdivide_indices(K.domain,h_BEM)
    N = length(Lₕ)
    M = length(Γ.IFS)
    #initialise Galerkin matrix:
    Galerkin_matrix = zeros(Complex{Float64},N,N)
    # create blank matrix of flags, describing if the matrix entry has been filled
    BEM_filled = zeros(Bool,N,N)
    m_count = 0


    # # save time computing repeated diagonal values if Γ is a uniform attactor
    # if Γ.uniform
    #     m = Lₕ[1]
    #     diag_val = singular_elliptic_double_integral(K, h_quad, m; Cosc=Cosc)
    # end

    if vary_quad
        h_quad_adjust = get_quad_scales(K::BIO,Lₕ::Array{Array{Int64,1},1})
    else
        h_quad_adjust = ones(length(Lₕ),length(Lₕ))
    end

    if Γ.uniform & repeat_blocks
        ℓ = length(Lₕ[1])
        diag_block_sizes = M.^(0:(ℓ-1))
    else
        diag_block_sizes = []
    end

    @showprogress 1 "Constructing BEM system " for m_count=1:length(Lₕ)#m in Lₕ
        # m_count+=1
        m = Lₕ[m_count]
        Γₘ = SubAttractor(Γ,m)
        # n_count = 0

        # if matrix is symmetric, will only need to compute ≈ half entries
        K.self_adjoint ? n_count_start = m_count : n_count_start = 1

        for n_count = n_count_start:length(Lₕ)#n in Lₕ
            if !BEM_filled[m_count,n_count] # check matrix entry hasn't been filled already
                # n_count += 1
                n = Lₕ[n_count]
                Γₙ = SubAttractor(Γ,n)
                x,y,w = barycentre_rule(Γₘ,Γₙ,h_quad*h_quad_adjust[m_count,n_count])
                if n==m
                    # if Γ.uniform
                    #     Galerkin_matrix[m_count,n_count] = diag_val
                    # else
                        Galerkin_matrix[m_count,n_count] = singular_elliptic_double_integral(K,h_quad,n;Cosc=Cosc)
                    # end
                else
                    Galerkin_matrix[m_count,n_count] = sum(w.*K.kernel(x,y))
                end

                # if matrix is symmetric, expoit this to save time
                if K.self_adjoint && n!=m
                    Galerkin_matrix[n_count,m_count] = Galerkin_matrix[m_count,n_count]
                end
            end
        end

        # now repeat entries along block diagonal, if at the end of a diagonal block, and uniform
        if in(m_count,diag_block_sizes)
            block_power = indexin(m_count,diag_block_sizes)[1]
            block_size = M^(block_power-1)
            for j=2:M
                block_range = ((j-1)*block_size+1):(j*block_size)
                Galerkin_matrix[block_range,block_range] = Galerkin_matrix[1:block_size,1:block_size]
                BEM_filled[block_range,block_range] .= 1
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
    return Projection(K.BIO.domain, K.Lₕ, coeffs)
end

struct Projection # onto the coefficient space of piecewise constants on the fractal
    domain::Fractal
    Lₕ::Array{Array{Int64,1},1} # subindices list
    coeffs::Array{<:Complex{<:Real},1}
end

"""
    SingleLayer(Γ::Fractal, k::Real=0.0)

represents the single layer boundary integral operator, Sϕ(x) = ∫_Γ Φ(x,y) ϕ(x) dHᵈ(y),
where Φ is the fundamental solution for the underlying PDE.
"""
function SingleLayer(Γ::Fractal, k::Number=0.0)
    if Γ.topological_dimension == 1
        if k==0.0 #2D Laplace case
            K = BIO(Γ, #fractal domain
            (x,y)->Φₜ(0.0,x,y), # log kernel
            (x,y)->zero_kernel(x,y), # kernel minus singularity
            0.0, # strength of singularity, corresponding to log singularity
            -1/(2π), # scaling of singularity
            true, #self-adjoint
            k #wavenumber
            )
        else #2D Helmholtz case        
            K = BIO(Γ, #fractal domain
            (x,y)->HelhmoltzGreen2D(k,x,y), # Hankel function
            (x,y)->HelhmoltzGreen2D_Lipschitz_part(k,x,y), # kernel minus singularity
            0.0, # strength of singularity, corresponding to log singularity
            -1/(2π), # scaling of singularity
            true, #self-adjoint
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
            true, #self-adjoint
            k #wavenumber
            )
        else #3D Helmholtz case        
            K = BIO(Γ, #fractal domain
            (x,y)->HelhmoltzGreen3D(k,x,y), # Green's function
            (x,y)->HelhmoltzGreen3D_Lipschitz_part(k,x,y), # kernel minus singularity
            1.0, # strength of singularity, corresponding to 1/|x-y|
            1/(4π), # scaling of singularity
            true, #self-adjoint
            k #wavenumber
            )
        end
    else
        error("Haven't coded single layer BIO for this many dimensions")
    end
end

function singular_elliptic_double_integral(K::BIO,h_quad::Real,index::Array{Int64}=[0]; Cosc = 2π, data=false)
    if real(K.wavenumber)*K.domain.diameter > Cosc
        return singular_elliptic_double_integral_full(K::BIO,h_quad::Real,index; Cosc = Cosc, data=data)
    else
        return singular_elliptic_double_integral_basic(K::BIO,h_quad::Real,index; data=data)
    end
end

function singular_elliptic_double_integral(Γ::Attractor,k::Number,h_quad::Real; Cosc = 2π, data=false)
    K = SingleLayer(Γ, k)
    singular_elliptic_double_integral(K, h_quad; Cosc = Cosc, data=data)
end

function singular_elliptic_double_integral_full(K::BIO,h_quad::Real,index::Array{Int64}=[0]; Cosc = 2π, data=false)
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
                if data
                    I_, Npts_ = eval_green_double_integral(Γₘ,K.singularity_strength,h; data=true)
                    I += K.singularity_scale*I_
                    Npts+=Npts_
                else
                    I += K.singularity_scale*eval_green_double_integral(Γₘ,K.singularity_strength,h; data=false)
                end
                I += sum(w.*K.Lipschitz_part_of_kernel(x,y))
            else
                I += sum(w.*K.kernel(x,y))
            end
        end
    end
    if data
        return I, Npts
    else
        return I
    end
end

function singular_elliptic_double_integral_basic(K::BIO,h::Real,index::Array{Int64}=[0]; data=false)
    Γ = SubAttractor(K.domain,index)
    x,y,w = barycentre_rule(Γ,Γ,h)
    if data
        Npts = length(w)
        I_, Npts_ = eval_green_double_integral(Γ,K.singularity_strength,h; data=true)
        Npts += Npts_
        I = K.singularity_scale*I_ + sum(w.*K.Lipschitz_part_of_kernel(x,y))
        return I, Npts
    else
        I = K.singularity_scale*eval_green_double_integral(Γ,K.singularity_strength,h) + sum(w.*K.Lipschitz_part_of_kernel(x,y))
        return I
    end
end

function single_layer_potential(Γ::Fractal,k::Real,density::Projection,points::Array{<:Real,2};h=0.1/k)
    num_points = size(points)[1]
    vals = zeros(Complex{Float64},num_points)
    # kernel is the same as the single layer BIO, so use that:
    S = SingleLayer(Γ, k)
    Φ = S.kernel
    for m_count = 1:length(density.Lₕ)
        m = density.Lₕ[m_count]
        Γₘ = SubAttractor(Γ,m)
        y,w = barycentre_rule(Γₘ,h)
        # now convert y to n+1 dimensions
        Y = [y zeros(length(w))]
        for n = 1:num_points
            vals[n] += w'*Φ(Y,points[n,:])*density.coeffs[m_count]
        end
    end
    return vals
end
# simplify for Laplace kernel case:
single_layer_potential(Γ::Fractal,density::Projection,points::Array{<:Real,2}) = single_layer_potential(Γ,0.0,density,points)

#have not tested far-field pattern yet:
function far_field_pattern(Γ::Fractal,k::Real,density::Projection,points::Array{<:Real,2};h=0.1/k)
    num_points = size(points)[1]
    vals = zeros(Complex{Float64},num_points)
    # far-field kernel:
    K(θ::Float64, ψ::Float64, x_1::Array{<:Real,1}, x_2::Array{<:Real,1}) = exp.(-im*k*(cos(θ)*x_1.+sin(ψ)*x_2))
    for m_count =1:length(density.Lₕ)
        m = density.Lₕ[m_count]
        Γₘ = SubAttractor(Γ,m)
        y,w = barycentre_rule(Γₘ,h)
        # now convert y to n+1 dimensions
        for n = 1:num_points
            vals[n] += w'*K(points[n,1],points[n,2],y[:,1],y[:,2])*density.coeffs[m_count]
        end
    end
    return vals
end

function far_field_pattern(Γ::Fractal,k::Real,density::Projection,points::Array{<:Real,1};h=0.1/k)
    num_points = length(points)
    vals = zeros(Complex{Float64},num_points)
    # far-field kernel:
    K(θ::Float64, x::Array{Float64,1}) = exp.(-im*k*(cos(θ)*x))
    for m_count =1:length(density.Lₕ)
        m = density.Lₕ[m_count]
        Γₘ = SubAttractor(Γ,m)
        y,w = barycentre_rule(Γₘ,h)
        # now convert y to n+1 dimensions
        for n = 1:num_points
            vals[n] += w'*K(points[n],y)*density.coeffs[m_count]
        end
    end
    return vals
end

# now some functions related to projections, and embeddings

function embed(f,g)
    #    make sure f is embedded into g
    if length(f.coeffs)>length(g.coeffs)
        F=f
        G=g
        f=G
        g=F
    end
    
    f_dim = length(f.coeffs)
    g_dim = length(g.coeffs)
    new_coeffs = zeros(Complex{Float64},g_dim)
    m_count = 1
    for m in g.Lₕ
        # find the vector in f
        n_count = 0
        for n in f.Lₕ
            n_count +=1
            if n == m[1:length(n)]
                new_coeffs[m_count] = f.coeffs[n_count]
                m_count+=1
                break
            end
        end
    end
    return Projection(f.domain,g.Lₕ,new_coeffs)
end

function -(f::Projection,g::Projection)
    ϕ = embed(f,g)
    if length(f.coeffs)>length(g.coeffs)
        return Projection(ϕ.domain,ϕ.Lₕ,f.coeffs-ϕ.coeffs)
    else
        return Projection(ϕ.domain,ϕ.Lₕ,ϕ.coeffs-g.coeffs)
    end
end

# The next couple of functions are designed to use less quadrature points in the BEM elements
# where this won't affect the accuracy.

F_nomeasure(r::Real, k::Number, n::Int64) = (1+(abs(k)*r)^(n/2+1))/r^(n+1)

function get_quad_scales(K::BIO,Lₕ::Array{Array{Int64,1},1})
    Γ = K.domain
    # compute upper and lower bounds for the F in my notes, which is stated above.
    F_upper = ones(length(Lₕ),length(Lₕ))
    F_lower = ones(length(Lₕ),length(Lₕ))
    for m_count = 1:length(Lₕ)
        Γₘ = SubAttractor(Γ,Lₕ[m_count])
        for n_count = m_count:length(Lₕ)
            Γₙ = SubAttractor(Γ,Lₕ[n_count])
            if n_count!=m_count
                dist_upper = dist(Γₙ.barycentre,Γₘ.barycentre)
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
