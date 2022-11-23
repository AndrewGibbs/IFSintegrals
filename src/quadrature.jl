function get_quadrature_from_partition(γ::partition_data, M::Int64, S::Vector{Similarity{T,M_}}, weights::Vector{Float64}, h::Float64) where {T<:Real, M_<:Real}
    if h<γ.diameter
        X = create_partition(γ, M, S, weights, h)
        N = length(X)
        w = zeros(N)
        zero_node = zero(T)
        x = [zero_node for _=1:N]
        for n=1:N
            x[n] = X[n].barycentre
            w[n] = X[n].weight
        end
    else
        x = [γ.barycentre]
        w = [γ.weight]
    end
    return x, w
end

function get_quadrature_from_partition(γ::partition_data, M::Int64, S::Vector{Similarity{T,M_}}, weights::Vector{Float64}, h::Float64) where {T<:AbstractVector, M_<:Union{Real,AbstractMatrix}}
    if h<γ.diameter
        X = create_partition(γ, M, S, weights, h)
        N = length(X)
        w = zeros(N)
        zero_node = zeros(Float64,length(γ.barycentre))
        x = [zero_node for _=1:N]
        for n=1:N
            x[n] = Vector{Float64}(X[n].barycentre)
            w[n] = X[n].weight
        end
    else
        x = [Vector{Float64}(γ.barycentre)]
        w = [γ.weight]
    end
    return x, w
end

"""
    x,w = barycentre_rule(Γ::Union{InvariantMeasure,SubInvariantMeasure},h::Real) 

returns a vector of N weights wⱼ>0 and nodes xⱼ ∈ Rⁿ, for approximation of integrals defined on an IFS Γ⊂Rⁿ.
"""
function barycentre_rule(Γ::InvariantMeasure{T,M_},h::Float64) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
    if Γ.homogeneous
        return get_quadrature_from_partition(partition_data_unindexed{T}(Γ.barycentre,Γ.measure,Γ.diameter), length(Γ.IFS), Γ.IFS, Γ.weights, h)
    else
        return get_quadrature_from_partition(partition_data_with_IFS{T,M_}(Γ.barycentre, Γ.measure, Γ.diameter, Γ.IFS), length(Γ.IFS), Γ.IFS, Γ.weights, h)
    end
end

function barycentre_rule(Γ::SubInvariantMeasure{T,M_},h::Float64) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}} 
    if Γ.parent_measure.homogeneous
        return get_quadrature_from_partition(partition_data_unindexed{T}(Γ.barycentre,Γ.measure,Γ.diameter), length(Γ.IFS), Γ.IFS, Γ.parent_measure.weights, h)
    else
        return get_quadrature_from_partition(partition_data_with_IFS{T,M_}(Γ.barycentre,Γ.measure,Γ.diameter,Γ.IFS), length(Γ.IFS), Γ.IFS, Γ.parent_measure.weights, h)
    end
end

"""
    x,y,w = barycentre_rule(Γ₁::Union{InvariantMeasure,SubInvariantMeasure},Γ₂::Union{InvariantMeasure,SubInvariantMeasure},h::Real)

returns N weights wⱼ>0 and nodes x,y ∈ Rⁿ, for approximation of double integrals over Γ₁,Γ₂⊂Rⁿ.
Uses Barycentre rule quadrature, the fractal Γ will be subdivided until each subcomponent has a diameter of
less than h.
"""
function barycentre_rule(Γ1::SelfSimilarFractal,Γ2::SelfSimilarFractal,h::Float64)
        x1, w1 = barycentre_rule(Γ1,h)
        n1 = length(w1)
        x2, w2 = barycentre_rule(Γ2,h)
        n2 = length(w2)
    return repeat(x1,inner=n2), repeat(x2,outer=n1), repeat(w1,inner=n2).*repeat(w2,outer=n1)
end

"""
    long_bary(Γ::SelfSimilarFractal{V,U},f::Function,h::Float64) where {V<:Union{Real,AbstractVector}, U<:Union{Real,AbstractMatrix}}

Computes the integral ∫_Γ f(x) dx via a barycentre rule, without storing all of the points.
This routine is slower than barycentre_rule, but is capable of achieving higher accuracy.
"""
function long_bary(Γ::SelfSimilarFractal{V,U},f::Function,h::Float64) where {V<:Union{Real,AbstractVector}, U<:Union{Real,AbstractMatrix}}
    # note that the input Γ may be an parent_measure or a subcomponent of an parent_measure
    I = 0.0 # value of integral which will be added to cumulatively
    N = 0
    M = length(Γ.IFS)
    if Γ.diameter>h
        for m=1:M
            I_,N_ = long_bary(Γ[m],f,h) # Here Γ[m] represents the m'th subcomponent of Γ, m=1,...,M#
            N += N_
            I += I_
        end
    else
        I = f(Γ.barycentre)*Γ.measure # one point quadrature at the finest level
        N = 1
    end
    return I, N
end

function long_bary(Γ₁::SelfSimilarFractal{V,U},Γ₂::SelfSimilarFractal{V,U},f::Function,h::Float64) where {V<:Union{Real,AbstractVector}, U<:Union{Real,AbstractMatrix}}
    # note that the input Γ may be an parent_measure or a subcomponent of an parent_measure
    I = 0.0 # value of integral which will be added to cumulatively
    N = 0
    M₁ = length(Γ₁.IFS)
    M₂ = length(Γ₂.IFS)
    if Γ₁.diameter>h && Γ₂.diameter>h
        for m=1:M₁
            for m_=1:M₂
                I_,N_ = long_bary(Γ₁[m],Γ₂[m_],f,h) # Here Γ[m] represents the m'th subcomponent of Γ, m=1,...,M#
                N += N_
                I += I_
            end
        end
    elseif Γ₁.diameter>h
        for m=1:M₁
            I_,N_ = long_bary(Γ₁[m],Γ₂,f,h) 
            N += N_
            I += I_
        end
    elseif Γ₂.diameter>h
        for m=1:M₂
            I_,N_ = long_bary(Γ₁,Γ₂[m],f,h)
            N += N_
            I += I_
        end
    else
        I = f(Γ₁.barycentre,Γ₂.barycentre)*Γ₁.measure*Γ₂.measure # one point quadrature at the finest level
        N = 1
    end
    return I, N
end

function long_bary_v2(Γ₁::SelfSimilarFractal{V,U},Γ₂::SelfSimilarFractal{V,U},f::Function,h::Float64) where {V<:Union{Real,AbstractVector}, U<:Union{Real,AbstractMatrix}}
    # note that the input Γ may be an parent_measure or a subcomponent of an parent_measure
    I = 0.0 # value of integral which will be added to cumulatively
    N = 0
    L₁ = subdivide_indices(Γ₁,h)
    L₂ = subdivide_indices(Γ₂,h)
    M₁ = [Γ₁[m] for m ∈ L₁]
    M₂ = [Γ₂[m] for m ∈ L₂]
    N₁ = length(M₁)
    N₂ = length(M₂)
    for n₁ = 1:N₁
        for n₂ =1:N₂
            I += f(M₁[n₁].barycentre, M₂[n₂].barycentre)*M₁[n₁].measure*M₂[n₂].measure
        end
    end
    return I
end

"""
    chaos_quad(Γ::SelfSimilarFractal{V,M},N::Int64;x₀=Γ.barycentre:::AbstractVector) where {V<:AbstractVector, M<:Union{Real,AbstractMatrix}}

Returns a vector of N weights w>0 and nodes x ∈ Rⁿ, for approximation of integrals defined on an IFS Γ⊂Rⁿ.
Using Chaos game quadrature.
Optional third input is the initial guess, which is taken as barycentre by default.
"""
function chaos_quad(Γ::SelfSimilarFractal{V,M}, N::Int64; x₀=Γ.barycentre::V) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    x = [Γ.barycentre for _ = 1:N]
    # x[1] = x₀
    for n=2:N
        x[n] = chaos_step(Γ,x[n-1])
    end
    return x, Γ.measure*ones(N)./N
end

function chaos_quad(Γ₁::SelfSimilarFractal{V,M},Γ₂::SelfSimilarFractal{V,M}, N::Int64) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    X₁,W₁ = chaos_quad(Γ₁,N)
    X₂,W₂ = chaos_quad(Γ₂,N)
    n1 = length(W₁)
    n2 = length(W₂)
    return repeat(X₁,inner=n2), repeat(X₂,outer=n1), repeat(W₁,inner=n2).*repeat(W₂,outer=n1)
end

function chaos_step(Γ::SelfSimilarFractal,x::Union{AbstractVector,Real})
    # randomly choose a map
    τ = minimum((1:length(getweights(Γ)))[rand().<cumsum(getweights(Γ))])
    return sim_map(Γ.IFS[τ],x)
end


"""
    x,w = gauss_quad(Γ::SelfSimilarFractal{V,M}, N::Int64) where {V<:Real, M<:Real}

Returns N Gaussian weights w ∈ Rᴺ and nodes x ∈ Rᴺˣᴺ.
Here Γ must be an SelfSimilarFractal in one spatial dimension.
N is the order of the Gauss rule, i.e. number of weights and nodes.
"""
function gauss_quad(Γ::InvariantMeasure{V,M}, N::Int64) where {V<:Real, M<:Real}
    J = get_Jacobi_matrix(Γ,N-1)
    vv = real.(eigvecs(J))
    x = real.(eigvals(J))
    w = vv[1,:].^2
    return x,w
end

# Bodged function below.
# I think there's a scaling issue when applying Mantica's algorithm to subcomonents.
 # A₀ = μ₁ does not seem to hold in that case. Need to understand why.
function gauss_quad(Γ::SubInvariantMeasure{V,M}, N::Int64) where {V<:Real, M<:Real}
    x, w = gauss_quad(Γ.parent_measure,N)
    sₘ = sim_comp(Γ.IFS, Γ.index)
    return [sim_map(sₘ,x[n]) for n=1:length(x)], w*Γ.measure/Γ.parent_measure.measure
end

function gauss_quad(Γ1::SelfSimilarFractal,Γ2::SelfSimilarFractal, N::Int64)
    x1, w1 = gauss_quad(Γ1, N)
    n1 = length(w1)
    x2, w2 = gauss_quad(Γ2, N)
    n2 = length(w2)
return repeat(x1,inner=n2), repeat(x2,outer=n1), repeat(w1,inner=n2).*repeat(w2,outer=n1)
end