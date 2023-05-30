abstract type DomainOperator{V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
end

"""
Define identity operator
"""
struct IdentityOperator{V,M} <: DomainOperator{V,M}# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    domain::SelfSimilarFractal{V,M}
    λ::Number
end

# scalar multiplication of identity
function *(c::Number, K::IdentityOperator{V,M}) where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    return IdentityOperator{V,M}(K.domain, c*K.λ)
end

# simplified constructor of identity operator
function Id(Γ::SelfSimilarFractal{V,M}) where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    return IdentityOperator(Γ,1.0)
end

function get_mass_matrix(I_Γ::IdentityOperator{V,M}, meshes::Vector{Vector{SubInvariantMeasure{V,M}}}) where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    # Γ = I_Γ.domain
    # Lₕ = subdivide_indices(Γ,h_mesh) #get vector indices for subcomponents
    N = sum([length(mesh) for mesh ∈ meshes])
    cum_mesh_size, total_num_els, _ = get_multi_mesh_sizes(meshes)
    # mesh = [SubInvariantMeasure(Γ,Lₕ[n]) for n=1:N] # partition Γ into subcomponents to make the mesh
    mass_matrix = zeros(total_num_els,total_num_els)
    for m=1:length(meshes)
        for n=1:cum_mesh_size[m]
            mass_matrix[cum_mesh_size[m]+n,cum_mesh_size[m]+n] = meshes[m][n].measure
        end
    end

    return I_Γ.λ * mass_matrix
end

"""
SIO is the type for singular integral operators.
"""
struct SIO{V,M} <: DomainOperator{V,M} where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    domain::Vector{SelfSimilarFractal{V,M}}
    kernel::Function
    Lipschitz_part_of_kernel::Function
    singularity_strength::Float64
    singularity_scale::Complex{<:Float64}
    self_adjoint::Bool
    wavenumber::Number
end

SIO(domain::Vector{SelfSimilarFractal{V,M}},
    kernel::Function,
    Lipschitz_part_of_kernel::Function,
    singularity_strength::Float64,
    singularity_scale::Complex{<:Float64},
    self_adjoint::Bool,
    wavenumber::Number
    ) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}} = SIO(
        [domain], kernel, Lipschitz_part_of_kernel, singularity_strength,
        singularity_scale, self_adjoint, wavenumber)

#constructor for zero wavenumber case
# SIO(domain::SelfSimilarFractal{V,M},kernel::Function,Lipschitz_part_of_kernel::Function,singularity_strength::Real,
# singularity_scale::Complex{<:Real},self_adjoint::Bool) where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}} =
# SIO{V,M}(domain,kernel,Lipschitz_part_of_kernel,singularity_strength,singularity_scale,self_adjoint,0.0)

"""
    DiscreteSIO(SIO::SIO; h_mesh::Real, h_quad::Real, h_quad_diag::Real)
    
is the constructor for a discretisation of a singular integral operator, 'SIO'.
h_mesh is the meshwidth parameter for the discretisation of the underlying fractal
h_quad denotes the discretisation parameter for the integrals in the stiffness matrix.
h_quad_diag is the parameter used to compute the diagonal elements of the matrix
"""
struct DiscreteSIO{V,M} <: DomainOperator{V,M} where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    domain::SelfSimilarFractal{V,M}
    SIO::DomainOperator{V,M}
    h_mesh::Float64
    h_quad::Float64
    mesh::Vector{Vector{SubInvariantMeasure{V,M}}} # one mesh per attractor
    Lₕ::Vector{Vector{Vector{Int64}}} # one set of subindices per attractor
    Galerkin_matrix::Matrix{<:Complex{<:Float64}} # eventually this type should be generalised
end

function check_for_similar_singular_integrals(Γ, prepared_singular_inds, n, n_)
    symmetry_group = get_symmetry_group(Γ)
    if Γ.disjoint
        # then we can immediately decide if integral is singular
        if n==n_
            is_similar_to_singular_integral = true
            n  != [0] ? ρ = 1/prod([Γ.IFS[nᵢ].r for nᵢ∈n]) : ρ = 1.0
            similar_index = 1
        else
            is_similar_to_singular_integral = false
            similar_index = Int64[]
            ρ = Float64[]
        end
    else # need to check if integral is similar to a singular integral, slightly longer process
        is_similar_to_singular_integral, ρ, similar_index = check_for_similar_integrals(Γ, prepared_singular_inds, n, n_, symmetry_group, symmetry_group, true)
    end
    return is_similar_to_singular_integral, ρ, similar_index
end

#constructor:#
function DiscreteSIO(K::SIO{V,M}; h_mesh::Real=max(2π/(10.0*K.wavenumber)), kwargs...) where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    Lₕ = subdivide_indices(K.domain,h_mesh) #get vector indices for subcomponents
    N = length(Lₕ)
    mesh = [SubInvariantMeasure(K.domain,Lₕ[n]) for n=1:N]
    return DiscreteSIO(K, mesh, Lₕ, h_mesh=h_mesh)
end

function DiscreteSIO(K::SIO{V,M_}, mesh::Vector{SubInvariantMeasure{V,M_}}, Lₕ::Vector{Vector{Int64}};
        h_mesh::Real= maximum([m.diameter for m∈mesh]),#max(2π/(10.0*K.wavenumber),K.domain.diameter+eps(),
        h_quad::Real=h_mesh, h_quad_diag::Real = h_quad,
        vary_quad::Bool = true, repeat_blocks::Bool =true,
        adjacency_function::Union{Function,Nothing}=nothing) where {V<:Union{Real,AbstractVector},M_<:Union{Real,AbstractMatrix}}
    Γ = K.domain
    # Lₕ = subdivide_indices(K.domain,h_mesh) #get vector indices for subcomponents
    N = length(Lₕ)
    # mesh = [SubInvariantMeasure(Γ,Lₕ[n]) for n=1:N] # partition Γ into subcomponents to make the mesh
    M = length(Γ.IFS)
    # symmetry_group = get_symmetry_group(Γ)
    # create blank matrix of flags, describing if the matrix entry has been filled
    BEM_filled = zeros(Bool,N,N)
    m_count = 0

    # Josh's stuff. Check if adjacency_function has been provided, if so, use it:
    if adjacency_function !== nothing
        use_users_adj_fn = true
        quad_type, quad_scale = adjacency_function(mesh)
    else
        use_users_adj_fn = false
    end


    # now get matrix of how much we can adjust the h_quad parameter for far away elements
    if vary_quad
        h_quad_adjust = get_quad_scales(K,Lₕ)
    else
        h_quad_adjust = ones(Float16, length(Lₕ),length(Lₕ))
    end

    if Γ.homogeneous && repeat_blocks && N>1
        ℓ = length(Lₕ[1])
        # the sizes of the repeated sub-blocks will be as follows:
        diag_block_sizes = M.^(0:(ℓ-1))
    else
         diag_block_sizes = []
    end

    #initialise Galerkin matrix:
    Galerkin_matrix = zeros(Complex{Float64},N,N)

    #collect the different types of singular matrix, and their indices
    h_high_scale = Γ.diameter*h_quad_diag/h_mesh
    s = K.singularity_strength
    # NOTE: when the non-disjoint singular stuff is optimised,
    # the first part of the below if statement can be replaced by the second
    if Γ.disjoint
        prepared_singular_inds = [([0],[0])]
        prepared_singular_vals = eval_green_double_integral(Γ, s, h_high_scale)
    else #non-disjoint
        A,B,prepared_singular_inds,R,log_stuff = construct_singularity_matrix(Γ, s)
        r = zeros(length(R))
        for n=1:length(r)
            (m,m_) = R[n]
            x,y,w = barycentre_rule(Γ[m],Γ[m_],h_high_scale)
            r[n] = w'*Φₜ.(s,x,y)
        end
        prepared_singular_vals = A\(B*r + log_stuff) # vector of 'singular values'
    end

    function scaler(ρ::Float64, m::Vector{<:Int64},m_::Vector{<:Int64},n::Vector{<:Int64},n_::Vector{<:Int64})
        # account for convention Γ₀:=Γ
        m  != [0] ? pₘ = prod(Γ.weights[m]) : pₘ = 1.0
        m_ != [0] ? pₘ_ = prod(Γ.weights[m_]) : pₘ_ = 1.0
        n  != [0] ? pₙ = prod(Γ.weights[n]) : pₙ = 1.0
        n_ != [0] ? pₙ_ = prod(Γ.weights[n_]) : pₙ_ = 1.0
        return ρ^(-s)*pₘ*pₘ_/pₙ/pₙ_
    end

    @showprogress 1 "Constructing discrete system " for m_count=1:N#m in Lₕ
        m = Lₕ[m_count]
        Γₘ = mesh[m_count]

        # if matrix is symmetric, will only need to compute ≈ half entries
        K.self_adjoint ? n_count_start = m_count : n_count_start = 1

        for n_count = n_count_start:N#n in Lₕ
            if !BEM_filled[m_count,n_count] # check matrix entry hasn't been filled already
                n = Lₕ[n_count]
                Γₙ = mesh[n_count] # mesh element

                # TO DO: add some disjointness condition here
                if use_users_adj_fn
                    similar_index = quad_type[n_count,m_count]
                    similar_index>0 ? is_similar_to_singular_integral = true : is_similar_to_singular_integral = false
                    ρ = quad_scale[n_count, m_count]
                    _, ρ_, _ = check_for_similar_singular_integrals(Γ, prepared_singular_inds, n, m)
                else
                    is_similar_to_singular_integral, ρ, similar_index = check_for_similar_singular_integrals(Γ, prepared_singular_inds, n, m)
                end

                if is_similar_to_singular_integral
                    similar_indices = prepared_singular_inds[similar_index]
                    scale_adjust = 1/scaler(ρ, similar_indices[1], similar_indices[2], n, m)
                    x,y,w = barycentre_rule(Γₘ,Γₙ,h_quad)
                    Galerkin_matrix[m_count,n_count] = K.singularity_scale*prepared_singular_vals[similar_index]*scale_adjust+ w'*K.Lipschitz_part_of_kernel.(x,y)
                    if K.singularity_strength == 0
                        m  != [0] ? pₘ = prod(Γ.weights[m]) : pₘ = 1.0
                        n  != [0] ? pₙ = prod(Γ.weights[n]) : pₙ = 1.0
                        Galerkin_matrix[m_count,n_count] += K.singularity_scale*Γ.measure^2*log(1/ρ)*pₙ*pₘ # log constant adjustment
                   end
                else
                    x,y,w = barycentre_rule(Γₘ,Γₙ,h_quad*h_quad_adjust[m_count,n_count]) # get quadrature
                    Galerkin_matrix[m_count,n_count] = w'*K.kernel.(x,y) # evaluate non-diagonal Galerkin integral
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
    DiscreteSIO(Γ, K, h_mesh, h_quad, mesh, Lₕ, Galerkin_matrix)
end

function get_multi_mesh_sizes(meshes::Vector{Vector{SubInvariantMeasure{V,M_}}})
    cum_mesh_size = vcat(0,cumsum([length(mesh) for mesh ∈ meshes]))
    num_meshes = length(meshes)
    mesh_el_indices = [(cum_mesh_size[n]+1):cum_mesh_size[n+1] for n=1:num_meshes]
    total_num_els = mesh_el_indices[end][end]
    return cum_mesh_size, total_num_els, mesh_el_indices
end
# cum_mesh_size, total_num_els, mesh_el_indices = get_multi_mesh_sizes(meshes)

# multiple scattering
function DiscreteSIO(K::SIO{V,M_}, meshes::Vector{Vector{SubInvariantMeasure{V,M_}}}, Lₕs::Vector{Vector{Vector{Int64}}};
    h_mesh::Real= maximum([m.diameter for m∈mesh]),#max(2π/(10.0*K.wavenumber),K.domain.diameter+eps(),
    h_quad::Real=h_mesh, h_quad_diag::Real = h_quad,
    vary_quad::Bool = true, repeat_blocks::Bool =true,
    adjacency_function::Union{Function,Nothing}=nothing) where {V<:Union{Real,AbstractVector},M_<:Union{Real,AbstractMatrix}}
    
    # get mesh indices and related numbers
    cum_mesh_size, total_num_els, mesh_el_indices = get_multi_mesh_sizes(meshes)

    Galerkin_matrix = zeros(ComplexF64,total_num_els,total_num_els)
    
    for (i,mesh_i) ∈ enumerate(meshes)
        for (j,mesh_j) ∈ enumerate(meshes)
            if i==j # diagonal entry
                diag_DSIO = DiscreteSIO(K, mesh_i, Lₕs[i])
                Galerkin_matrix[mesh_el_indices[i], mesh_el_indices[j]] = diag_DSIO.Galerkin_matrix
            else # compute from scratch, simple quadrature
                for (𝙈ₘᵢ, m_count) ∈ enumerate(mesh_i)
                    for (𝙈ₙⱼ, n_count) ∈ enumerate(mesh_j)
                        x,y,w = barycentre_rule(𝙈ₘᵢ, 𝙈ₙⱼ, h_quad) # get quadrature
                        Galerkin_matrix[cum_mesh_size[i]+m_count, cum_mesh_size[j]+n_count] = w'*K.kernel.(x,y) # evaluate non-diagonal Galerkin integral
                    end
                end
            end
        end
    end
    DiscreteSIO(Γ, K, h_mesh, h_quad, mesh, Lₕ, Galerkin_matrix)
end

"""
    SingleLayer(Γ::SelfSimilarFractal, wavenumber::Real=0.0)

represents the single layer boundary integral operator, Sϕ(x) = ∫_Γ Φ(x,y) ϕ(x) dHᵈ(y),
where Φ is the fundamental solution for the underlying PDE.
"""
# function SingleLayer(Γ::SelfSimilarFractal{V,M}, k::Number=0.0) where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
#     if Γ.spatial_dimension == 1
#         if k==0.0 #2D Laplace case
#             K = SIO{V,M}(Γ, #fractal domain
#             (x,y)->Φₜ(0.0,x,y), # log kernel
#             (x,y)->zero_kernel(x,y), # kernel minus singularity
#             0.0, # strength of singularity, corresponding to log singularity
#             -1/(2π), # scaling of singularity
#             true, #self-adjoint
#             k #wavenumber
#             )
#         else #2D Helmholtz case        
#             K = SIO{V,M}(Γ, #fractal domain
#             (x,y)->HelhmoltzGreen2D(k,x,y), # Hankel function
#             (x,y)->HelhmoltzGreen2D_Lipschitz_part(k,x,y), # kernel minus singularity
#             0.0, # strength of singularity, corresponding to log singularity
#             -1/(2π), # scaling of singularity
#             true, #self-adjoint
#             k #wavenumber
#             )
#         end
#     elseif Γ.spatial_dimension == 2
#         if k==0.0 #3D Laplace case
#             K = SIO{V,M}(Γ, #fractal domain
#             (x,y)-> Φₜ(1.0,x,y), # Green's function
#             (x,y)-> zero_kernel(x,y), # kernel minus singularity
#             1.0, # strength of singularity, corresponding to 1/|x-y|
#             1/(4π), # scaling of singularity
#             true, #self-adjoint
#             k #wavenumber
#             )
#         else #3D Helmholtz case        
#             K = SIO{V,M}(Γ, #fractal domain
#             (x,y)->HelhmoltzGreen3D(k,x,y), # Green's function
#             (x,y)->HelhmoltzGreen3D_Lipschitz_part(k,x,y), # kernel minus singularity
#             1.0, # strength of singularity, corresponding to 1/|x-y|
#             1/(4π), # scaling of singularity
#             true, #self-adjoint
#             k #wavenumber
#             )
#         end
#     else
#         error("Haven't coded single layer SIO for this many dimensions")
#     end
# end

# function VolumePotential(Γ::SelfSimilarFractal{V,M}, k::Number) where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
#     if Γ.spatial_dimension == 2       
#         K = SIO{V,M}(Γ, #fractal domain
#         (x,y)->HelhmoltzGreen2D(k,x,y), # Hankel function
#         (x,y)->HelhmoltzGreen2D_Lipschitz_part(k,x,y), # kernel minus singularity
#         0.0, # strength of singularity, corresponding to log singularity
#         -1/(2π), # scaling of singularity
#         true, #self-adjoint
#         k #wavenumber
#         )
#     elseif Γ.spatial_dimension == 3    
#         K = SIO{V,M}(Γ, #fractal domain
#         (x,y)->HelhmoltzGreen3D(k,x,y), # Green's function
#         (x,y)->HelhmoltzGreen3D_Lipschitz_part(k,x,y), # kernel minus singularity
#         1.0, # strength of singularity, corresponding to 1/|x-y|
#         1/(4π), # scaling of singularity
#         true, #self-adjoint
#         k #wavenumber
#         )
#     else
#         error("Haven't coded single layer SIO for this many dimensions")
#     end
# end

struct OperatorSum{V,M} <: DomainOperator{V,M} where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    domain::Vector{SelfSimilarFractal{V,M}}
    operators::Vector{<:DomainOperator{V,M}}
end

function *(c::Number, K::SIO{V,M})# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    return SIO{V,M}(K.domain,
    (x,y) -> c*K.kernel(x,y),
    (x,y) -> c*K.Lipschitz_part_of_kernel(x,y),
    K.singularity_strength,
    c*K.singularity_scale,
    K.self_adjoint,
    K.wavenumber
    )
end

*(K::SIO, c::Number) = c*K
*(c::Number, K::OperatorSum) = OperatorSum(K.domain, [c*J for J∈K.operators])
*(K::OperatorSum, c::Number) = c*K
-(G::DomainOperator,F::DomainOperator) = G + (-1.0*F)

function +(F::DomainOperator{V,M}, G::DomainOperator{V,M}) where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    F.domain == G.domain ? nothing : error("domains of operators must match")
    return OperatorSum(F.domain, [F,G])
end

function +(F::DomainOperator{V,M}, G::OperatorSum{V,M}) where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    F.domain == G.domain ? nothing : error("domains of operators must match")
    return OperatorSum(F.domain, vcat([F],G.operators))
end

function +(G::OperatorSum{V,M},F::DomainOperator{V,M}) where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    F.domain == G.domain ? nothing : error("domains of operators must match")
    return OperatorSum(F.domain, vcat(G.operators,[F]))
end

function DiscreteSIO(K::OperatorSum; h_mesh::Real=K.domain.diameter, kwargs...)

    Γ = K.domain

    Lₕs  = [subdivide_indices(Γ,h_mesh) for Γ ∈ K.domain]# get vector indices for subcomponents
    meshes = [[SubInvariantMeasure(Γ[j],Lₕs[j][n]) for n=1:N] for j=1:length(domain)]
    # N = sum([length(Lₕ) for Lₕ ∈ Lₕs])

    _, total_num_els, _ = get_multi_mesh_sizes(meshes)

    # construct stiffness_matrix by summing together matrices from different discrete IOs and mass matrices
    stiffness_matrix = zeros(ComplexF64,total_num_els,total_num_els)

    for J ∈ K.operators
        if isa(J,SIO)
            stiffness_matrix += DiscreteSIO(J, meshes, Lₕs; kwargs...).Galerkin_matrix
        elseif isa(J,IdentityOperator)
            stiffness_matrix += get_mass_matrix(J, meshes)
        end
    end

    if haskey(kwargs,"h_quad")
        h_quad = kwargs["h_quad"]
    else
        h_quad = h_mesh
    end

    return DiscreteSIO(Γ, K, h_mesh, h_quad, mesh, Lₕ, stiffness_matrix)
end

function singular_elliptic_double_integral(K::SIO, h_quad::Real, index::Array{Int64}=[0]; Cosc = 2π)
    # if there are many wavelengths, use the more complicated approximation, which has better rates for high frequencies
    if real(K.wavenumber)*SubInvariantMeasure(K.domain,index).diameter > Cosc
        return singular_elliptic_double_integral_full(K::SIO,h_quad::Real,index; Cosc = Cosc)
    else
        return singular_elliptic_double_integral_basic(K::SIO,h_quad::Real,index)
    end
end

function singular_elliptic_double_integral(Γ::Union{InvariantMeasure,SubInvariantMeasure},k::Number,h_quad::Real; Cosc = 2π)
    K = SingleLayer(Γ, k)
    return singular_elliptic_double_integral(K, h_quad; Cosc = Cosc)
end

function singular_elliptic_double_integral_full(K::SIO,h_quad::Real,index::Array{Int64}=[0]; Cosc = 2π)
    # following the procedure in the paper, to avoid errors in far-field.
    # notation agrees with quadrature paper here, rather than BEM paper.
    Γ = SubInvariantMeasure(K.domain,index)
    h_star = Cosc/real(K.wavenumber) # new quadrature parameter introduced for this routine
    h = min(h_star,h_quad)
    L_h_star = subdivide_indices(Γ,h_star)
    Npts = 0

    I = zero(Complex{Float64})
    for n in L_h_star
        Γₙ = SubInvariantMeasure(Γ,n)
        for m in L_h_star
            Γₘ = SubInvariantMeasure(Γ,m)
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

function singular_elliptic_double_integral_basic(K::SIO,h::Real,index::Array{Int64}=[0])
    Γ = SubInvariantMeasure(K.domain,index)
    x,y,w = barycentre_rule(Γ,Γ,h)
    if K.domain.disjoint
        I = K.singularity_scale*eval_green_double_integral(Γ,K.singularity_strength,h) + w'*K.Lipschitz_part_of_kernel.(x,y)
    else
        I = K.singularity_scale*s_energy(Γ,K.singularity_strength,h) + w'*K.Lipschitz_part_of_kernel.(x,y)
    end
    return I
end

# The next couple of functions are designed to use less quadrature points in the BEM elements
# where this won't affect the accuracy.

F_nomeasure(r::Real, k::Number, n::Int64) = (1+(abs(k)*r)^(n/2+1))/r^(n+1)

function get_quad_scales(K::SIO,Lₕ::Vector{Vector{Int64}})
    Γ = K.domain
    # compute upper and lower bounds for the F in my notes, which is stated above.
    F_upper = ones(Float64,length(Lₕ),length(Lₕ))
    F_lower = ones(Float64,length(Lₕ),length(Lₕ))
    for m_count = 1:length(Lₕ)
        Γₘ = SubInvariantMeasure(Γ,Lₕ[m_count])
        for n_count = m_count:length(Lₕ)
            Γₙ = SubInvariantMeasure(Γ,Lₕ[n_count])
            if n_count!=m_count
                dist_upper = norm(Γₙ.barycentre-Γₘ.barycentre)
                dist_lower = max(dist_upper - Γₙ.diameter - Γₘ.diameter,0)
                measure_weight = Γₙ.measure*Γₘ.measure
                # noting that F_nomeasure is monotonic decreasing in r, can bound as follows:
                if dist_lower>0
                    F_upper[m_count,n_count] = measure_weight*F_nomeasure(dist_lower, K.wavenumber, K.domain.spatial_dimension)
                else
                    F_upper[m_count,n_count]= Inf
                end
                 F_lower[m_count,n_count] = measure_weight*F_nomeasure(dist_upper, K.wavenumber, K.domain.spatial_dimension)
            else
                F_upper[m_count,n_count] = Inf
            end
            
            # now by symmetry
            F_upper[n_count,m_count] = F_upper[m_count,n_count]
            F_lower[n_count,m_count] = F_lower[m_count,n_count]
        end
    end
    # now can get a lower bound estimate on the quantity from my notes:
    quad_scales = Float16.(floor.(max.(sqrt.(maximum(F_lower)./F_upper),1),sigdigits = 4))
    # have made a change which rounds down to nearest Float16, to save space
    return quad_scales
end
