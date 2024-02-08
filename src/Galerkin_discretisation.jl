"""
    DiscreteSIO(SIO::SIO; h_mesh::Real, h_quad::Real, h_quad_diag::Real)
    
is the constructor for a discretisation of a singular integral operator, 'SIO'.
h_mesh is the meshwidth parameter for the discretisation of the underlying fractal
h_quad denotes the discretisation parameter for the integrals in the stiffness matrix.
h_quad_diag is the parameter used to compute the diagonal elements of the matrix
"""
struct DiscreteSIO{V,M,Ω<:FractalMeasure{V,M},T<:Union{AbstractMatrix,LinearMap}}# <: DomainOperator{Ω<:FractalMeasure{V,M}}# where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    #domain::FractalMeasure{V,M}
    SIO::DomainOperator{Ω}
    h_mesh::Float64
    h_quad::Float64
    mesh::Vector{SubInvariantMeasure{V,M}} # one mesh per attractor
    Galerkin_matrix::T # eventually this type should be generalised
    mesh_el_indices::Vector{UnitRange{Int64}}
    use_gmres::Bool
end
# eventually the Galerkin_matrix type should be generalised, to <Union{AbstractMatrix,AbstractLinearOperator}

#(K, h_mesh, h_quad, meshes, Galerkin_matrix, mesh_el_indices)

function check_for_similar_singular_integrals(Γ::InvariantMeasure,
                                            prepared_singular_inds::Vector{Tuple{Vector{Int64}, Vector{Int64}}}, 
                                            n::Vector{Int64},
                                            n_::Vector{Int64})
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

function get_multi_mesh_sizes(meshes::Vector{Vector{SubInvariantMeasure{V,M_}}}
                        ) where {V<:Union{Real,AbstractVector},M_<:Union{Real,AbstractMatrix}}
    cum_mesh_size = vcat(0,cumsum([length(mesh) for mesh ∈ meshes]))
    num_meshes = length(meshes)
    mesh_el_indices = [(cum_mesh_size[n]+1):cum_mesh_size[n+1] for n=1:num_meshes]
    total_num_els = mesh_el_indices[end][end]
    return cum_mesh_size, total_num_els, mesh_el_indices
end
# cum_mesh_size, total_num_els, mesh_el_indices = get_multi_mesh_sizes(meshes)

# (single) singular operator on two disjoint scatterers
function discretise_Galerkin_block(
    K::SIO{Ω},
    mesh1::Vector{SubInvariantMeasure{V,M_}},
    mesh2::Vector{SubInvariantMeasure{V,M_}},
    h_quad::Number;
    kwargs...
    ) where 
    {V<:Union{Real,AbstractVector},
    M_<:Union{Real,AbstractMatrix},
    Ω<:UnionInvariantMeasure{V,M_}}

    # initialise Galerkin matrix
    Galerkin_matrix = zeros(ComplexF64,length(mesh1),length(mesh2))

    # get quad scales
    quad_scales = get_quad_scales(abs(K.wavenumber),
                    K.domain.spatial_dimension,
                    mesh1, mesh2)

    # construct Galerkin matrix
    @showprogress 1 "Constructing off-diagonal Galerkin block" for (m, 𝙈ₘ) ∈ enumerate(mesh1)
        for (n, Mₙ) ∈ enumerate(mesh2)
            x,y,w = barycentre_rule(𝙈ₘ, Mₙ, h_quad*quad_scales[m,n]) # get quadrature
            Galerkin_matrix[m, n] = w'*K.kernel.(x,y) # evaluate non-diagonal Galerkin integral
        end
    end
    return Galerkin_matrix
end

# sum of operators on two disjoint scatterers
function discretise_Galerkin_block(
    K::OperatorSum{Ω},
    mesh1::Vector{SubInvariantMeasure{V,M_}},
    mesh2::Vector{SubInvariantMeasure{V,M_}},
    h_quad::Number;
    kwargs...
    ) where 
    {V<:Union{Real,AbstractVector},
    M_<:Union{Real,AbstractMatrix},
    Ω<:FractalMeasure{V,M_}}

    # initialise Galerkin matrix
    Galerkin_matrix = zeros(ComplexF64,length(mesh1),length(mesh2))

    # now sum up all relevant Galerkin matrices for each operator
    for A ∈ K.operators
        if !(typeof(A)<:IdentityOperator) # identity operator will be zero
            Galerkin_matrix .+= 
                discretise_Galerkin_block(A, mesh1, mesh2, h_quad; kwargs...)
        end
    end

    return Galerkin_matrix
end

# sum of operators on single scatterer
function discretise_Galerkin_block(
    K::OperatorSum{Ω},
    mesh::Vector{SubInvariantMeasure{V,M_}},
    h_quad::Number;
    kwargs...
    ) where 
    {V<:Union{Real,AbstractVector},
    M_<:Union{Real,AbstractMatrix},
    Ω<:FractalMeasure{V,M_}}

    # initialise Galerkin matrix
    Galerkin_matrix = zeros(ComplexF64,length(mesh),length(mesh))

    # now sum up all relevant Galerkin matrices for each operator
    for A ∈ K.operators
        if isa(A,IdentityOperator) # identity operator will be zero
            Galerkin_matrix .+= get_mass_matrix(A,mesh)
        elseif isa(A,SIO)
            Galerkin_matrix .+= discretise_Galerkin_block(A,mesh,h_quad;kwargs...)
        end
    end

    return Galerkin_matrix
end

# single singular operator on single scatterer
function discretise_Galerkin_block(
    K::SIO{Ω},
    mesh::Vector{SubInvariantMeasure{V,M_}},
    h_quad::Number;
    h_quad_diag::Real = h_quad,
    vary_quad::Bool = true,
    repeat_blocks::Bool =true,
    adjacency_function::Union{Function,Nothing}=nothing#,
    # kwargs...
    ) where 
    {V<:Union{Real,AbstractVector},
    M_<:Union{Real,AbstractMatrix},
    Ω<:FractalMeasure{V,M_}}

    # focus on the fractal where the mesh lives
    Γ = mesh[1].parent_measure
    N = length(mesh)
    M = length(Γ.IFS)
    BEM_filled = zeros(Bool,N,N)
    m_count = 0

    # reconstruct h_mesh
    h_mesh = maximum([m.diameter for m ∈ mesh])

    # Josh's stuff. Check if adjacency_function has been provided, if so, use it:
    if adjacency_function !== nothing
        use_users_adj_fn = true
        quad_type, quad_scale = adjacency_function(mesh)
    else
        use_users_adj_fn = false
    end

    # now get matrix of how much we can adjust the h_quad parameter for far away elements
    if vary_quad
        h_quad_adjust = get_quad_scales(abs(K.wavenumber), Γ.spatial_dimension, mesh)
    else
        h_quad_adjust = ones(Float16, length(mesh),length(mesh))
    end

    if Γ.homogeneous && repeat_blocks && N>1
        ℓ = length(mesh[1].index)
        # the sizes of the repeated sub-blocks will be as follows:
        diag_block_sizes = M.^(0:(ℓ-1))
    else
         diag_block_sizes = []
    end

    #initialise Galerkin matrix:
    Galerkin_matrix = zeros(ComplexF64,N,N)

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

    @showprogress 1 "Constructing Galerkin block " for m_count=1:N#m in Lₕ
        m = mesh[m_count].index
        Γₘ = mesh[m_count]

        # if matrix is symmetric, will only need to compute ≈ half entries
        K.self_adjoint ? n_count_start = m_count : n_count_start = 1

        for n_count = n_count_start:N#n in Lₕ
            if !BEM_filled[m_count,n_count] # check matrix entry hasn't been filled already
                n = mesh[n_count].index
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
                        # println(Γ.measure^2*log(1/ρ)*pₙ*pₘ)
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
    return Galerkin_matrix
end

function mesh_fractal(Γ::InvariantMeasure, h_mesh::Number)
    Lₕ = subdivide_indices(Γ, h_mesh)
    return [ Γ[m] for m ∈ Lₕ]
end

# main entrance
function DiscreteSIO(
    K::DomainOperator{Ω};# meshes::Vector{Vector{SubInvariantMeasure{V,M_}}}, Lₕs::Vector{Vector{Vector{Int64}}};
    h_mesh::Number=Inf,
    h_quad = h_mesh,
    use_gmres = false,
    kwargs...
    ) where {V<:Union{Real,AbstractVector},
            M_<:Union{Real,AbstractMatrix},
            Ω<:FractalMeasure{V,M_}}
    
    if isa(K.domain,InvariantMeasure) # one mesh in collection of meshes
        meshes = [mesh_fractal(K.domain, h_mesh)]
        # block_message = "blocks completed"
    elseif isa(K.domain,UnionInvariantMeasure)
        meshes = [ mesh_fractal(Γ, h_mesh) for Γ ∈ K.domain.invariant_measures]
        # block_message = ""
    else
        error("domain type not supported, use InvariantMeasure instead")
    end
    
    # get mesh indices and related numbers
    cum_mesh_size, total_num_els, mesh_el_indices = get_multi_mesh_sizes(meshes)

    Galerkin_matrix = zeros(ComplexF64,total_num_els,total_num_els)
    
    for (i,mesh_i) ∈ enumerate(meshes)
        for (j,mesh_j) ∈ enumerate(meshes)
            if i==j # diagonal entry
                Galerkin_matrix[mesh_el_indices[i], mesh_el_indices[j]] = 
                    discretise_Galerkin_block(K, mesh_j, h_quad; kwargs...)
            elseif j<i && K.self_adjoint
                Galerkin_matrix[mesh_el_indices[i], mesh_el_indices[j]] = 
                    transpose(Galerkin_matrix[mesh_el_indices[j], mesh_el_indices[i]])
            else# j>i && K.self_adjoint# compute from scratch, simple quadrature
                Galerkin_matrix[mesh_el_indices[i], mesh_el_indices[j]] = 
                    discretise_Galerkin_block(K, mesh_i, mesh_j, h_quad; kwargs...)
            end
        end
    end
    DiscreteSIO{V,M_,Ω,Matrix{ComplexF64}}(K, h_mesh, h_quad, vcat(meshes...), Galerkin_matrix, mesh_el_indices, use_gmres)
end