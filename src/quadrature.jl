function get_quadrature_from_partition(γ::partition_data{T}, M::Int64, S::Vector{Similarity{T,M_}}, weights::Vector{Float64}, h::Float64) where {T<:Real, M_<:Real}
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

function get_quadrature_from_partition(γ::partition_data{T}, M::Int64, S::Vector{Similarity{T,M_}}, weights::Vector{Float64}, h::Float64) where {T<:AbstractVector, M_<:Union{Real,AbstractMatrix}}
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
    x,w = barycentre_rule(Γ::Union{Attractor,SubAttractor},h::Real) 

returns N weights w ∈ Rⁿ and nodes x ∈ Rᴺˣⁿ,
for approximation of integrals defined on an IFS Γ
"""
function barycentre_rule(Γ::Attractor{T,M_},h::Float64) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
    return get_quadrature_from_partition(partition_data_unindexed{T}(Γ.barycentre,Γ.measure,Γ.diameter), length(Γ.IFS), Γ.IFS, Γ.weights, h)
end

function barycentre_rule(Γ::SubAttractor{T,M_},h::Float64) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}} 
    return get_quadrature_from_partition(partition_data_unindexed{T}(Γ.barycentre,Γ.measure,Γ.diameter), length(Γ.IFS), Γ.IFS, Γ.attractor.weights, h)
end

"""
    x,y,w = barycentre_rule(Γ₁::Union{Attractor,SubAttractor},Γ₂::Union{Attractor,SubAttractor},h::Real)

returns N weights w ∈ Rⁿ and nodes x,y ∈ Rᴺˣⁿ,
for approximation of double integrals over Γ₁,Γ₂.
"""
function barycentre_rule(Γ1::SelfSimilarFractal,Γ2::SelfSimilarFractal,h::Float64)
        x1, w1 = barycentre_rule(Γ1,h)
        n1 = length(w1)
        x2, w2 = barycentre_rule(Γ2,h)
        n2 = length(w2)
        # X1 = [zero(x1[1]) for _=1:n1*n2]
        # X2 = [zero(x2[1]) for _=1:n1*n2]
        # W = zeros(Float64, n1*n2)
        # for i=0:(n1-1)
        #     for j=0:(n2-1)
        #         X1[j*n1 + i + 1] = x1[i + 1]
        #         X2[j*n1 + i + 1] = x2[j + 1]
        #         W[j*n1 + i + 1] = w1[i + 1]*w2[j + 1]
        #     end
        # end
    return repeat(x1,inner=n2), repeat(x2,outer=n1), repeat(w1,inner=n2).*repeat(w2,outer=n1)
end

function long_bary(Γ::SelfSimilarFractal{V,U},f::Function,h::Float64) where {V<:Union{Real,AbstractVector}, U<:Union{Real,AbstractMatrix}}
    # note that the input Γ may be an attractor or a subcomponent of an attractor
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

function chaos_quad(Γ::SelfSimilarFractal{V,M},N::Int64,x₀::AbstractVector) where {V<:AbstractVector, M<:Union{Real,AbstractMatrix}}
    x = [Vector{Float64}(undef, length(x₀)) for _ = 1:N]
    x[1] = x₀
    for n=2:N
        x[n] = chaos_step(Γ,x[n-1])
    end
    return x, ones(N)./N
end

function chaos_quad(Γ::SelfSimilarFractal{V,M},N::Int64,x₀::Float64) where {V<:Real, M<:Real}
    x = zeros(V,N)
    x[1] = x₀
    for n=2:N
        x[n] = chaos_step(Γ,x[n-1])
    end
    return x, ones(N)./N
end

function chaos_step(Γ::SelfSimilarFractal,x::Union{AbstractVector,Real})
    # randomly choose a map
    τ = minimum((1:length(Γ.weights))[rand().<cumsum(Γ.weights)])
    return sim_map(Γ.IFS[τ],x)
end