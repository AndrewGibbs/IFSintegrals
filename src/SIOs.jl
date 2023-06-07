abstract type DomainOperator{Ω<:FractalMeasure}# where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
end
# changing stuff
"""
Define identity operator
"""
struct IdentityOperator{Ω} <: DomainOperator{Ω}# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    domain::Ω#FractalMeasure{V,M}
    λ::Number
    self_adjoint::Bool
end

# scalar multiplication of identity
function *(c::Number, K::IdentityOperator)# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    return IdentityOperator(K.domain, c*K.λ, true)
end

# simplified constructor of identity operator
function Id(Γ::FractalMeasure)# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    return IdentityOperator(Γ,1.0, true)
end

function get_mass_matrix(I_Γ::IdentityOperator, mesh::Vector{<:SubInvariantMeasure})# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    N = length(mesh)
    # mesh = [SubInvariantMeasure(Γ,Lₕ[n]) for n=1:N] # partition Γ into subcomponents to make the mesh
    mass_matrix = zeros(N,N)
    for n=1:N
        mass_matrix[n,n] = mesh[n].measure
    end

    return I_Γ.λ * mass_matrix
end

"""
SIO is the type for singular integral operators.
"""
struct SIO{Ω} <: DomainOperator{Ω}# where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    domain::Ω#Vector{FractalMeasure{V,M}}
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
    self_adjoint::Bool
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
*(c::Number, K::OperatorSum) = OperatorSum(K.domain, [c*J for J∈K.operators], K.self_adjoint)
*(K::OperatorSum, c::Number) = c*K
-(G::DomainOperator,F::DomainOperator) = G + (-1.0*F)

function +(F::DomainOperator, G::DomainOperator)# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    F.domain == G.domain ? nothing : error("domains of operators must match")
    return OperatorSum(F.domain, [F,G], prod([F.self_adjoint,G.self_adjoint]))
end

function +(F::DomainOperator, G::OperatorSum)# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    F.domain == G.domain ? nothing : error("domains of operators must match")
    return OperatorSum(F.domain, vcat([F],G.operators), prod([F.self_adjoint,G.self_adjoint]))
end

function +(G::OperatorSum,F::DomainOperator)# where {V<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
    F.domain == G.domain ? nothing : error("domains of operators must match")
    return OperatorSum(F.domain, vcat(G.operators,[F]), prod([F.self_adjoint,G.self_adjoint]))
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