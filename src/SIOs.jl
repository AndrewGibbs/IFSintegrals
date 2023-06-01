abstract type DomainOperator{Ω<:SelfSimilarFractal}# where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
end
# changing stuff
"""
Define identity operator
"""
struct IdentityOperator{Ω} <: DomainOperator{Ω}# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    domain::Ω#SelfSimilarFractal{V,M}
    λ::Number
end

# scalar multiplication of identity
function *(c::Number, K::IdentityOperator)# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    return IdentityOperator(K.domain, c*K.λ)
end

# simplified constructor of identity operator
function Id(Γ::SelfSimilarFractal)# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    return IdentityOperator(Γ,1.0)
end

function get_mass_matrix(I_Γ::IdentityOperator, meshes::Vector{Vector{SubInvariantMeasure}})# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
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
struct SIO{Ω} <: DomainOperator{Ω}# where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    domain::Ω#Vector{SelfSimilarFractal{V,M}}
    kernel::Function
    Lipschitz_part_of_kernel::Function
    singularity_strength::Float64
    singularity_scale::ComplexF64
    self_adjoint::Bool
    wavenumber::Float64
end

struct OperatorSum{Ω} <: DomainOperator{Ω}# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    domain::Ω
    operators::Vector{<:DomainOperator{Ω}}
end

function *(c::Number, K::SIO)# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    return SIO(K.domain,
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

function +(F::DomainOperator, G::DomainOperator)# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    F.domain == G.domain ? nothing : error("domains of operators must match")
    return OperatorSum(F.domain, [F,G])
end

function +(F::DomainOperator, G::OperatorSum)# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    F.domain == G.domain ? nothing : error("domains of operators must match")
    return OperatorSum(F.domain, vcat([F],G.operators))
end

function +(G::OperatorSum,F::DomainOperator)# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    F.domain == G.domain ? nothing : error("domains of operators must match")
    return OperatorSum(F.domain, vcat(G.operators,[F]))
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

function get_quad_scales(k::Real, spat_dim::Int64, mesh::Vector{<:SubInvariantMeasure})
    # compute upper and lower bounds for the F in my notes, which is stated above.
    N = length(mesh)
    F_upper = ones(Float64,N,N)
    F_lower = ones(Float64,N,N)
    for m_count = 1:N
        Γₘ = mesh[m_count]
        for n_count = m_count:N
            Γₙ = mesh[n_count]
            if n_count!=m_count
                dist_upper = norm(Γₙ.barycentre-Γₘ.barycentre)
                dist_lower = max(dist_upper - Γₙ.diameter - Γₘ.diameter,0)
                measure_weight = Γₙ.measure*Γₘ.measure
                # noting that F_nomeasure is monotonic decreasing in r, can bound as follows:
                if dist_lower>0
                    F_upper[m_count,n_count] = measure_weight*F_nomeasure(dist_lower, k, spat_dim)
                else
                    F_upper[m_count,n_count]= Inf
                end
                 F_lower[m_count,n_count] = measure_weight*F_nomeasure(dist_upper, k, spat_dim)
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
