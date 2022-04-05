import Base: zero

abstract type partition_data{T<:Union{Real,AbstractVector}}
end

struct partition_data_indexed{T<:Union{Real,AbstractVector}} <: partition_data{T}
    barycentre::T
    weight::Float64
    diameter::Float64
    index::Vector{UInt64}
end

struct partition_data_unindexed{T<:Union{Real,AbstractVector}} <: partition_data{T}
    barycentre::T
    weight::Float64
    diameter::Float64
end

zero(::Type{partition_data_indexed{T}}) where {T<:Union{Real,AbstractVector}} =  partition_data_indexed{T}(zero(T),0.0,0.0,[0])

zero(::Type{partition_data_unindexed{T}}) where {T<:Union{Real,AbstractVector}} =  partition_data_unindexed{T}(zero(T),0.0,0.0)

function subdivide_indexed(M::Integer,S::Vector{Similarity{T,M_}},weights::Vector{Float64}, X::partition_data_indexed{T}) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
    Y = zeros(partition_data_indexed{T},M)
    for (m,sₘ) in enumerate(S)
        Y[m] = partition_data_indexed{T}(sim_map(sₘ,X.barycentre), X.weight*weights[m], X.diameter*sₘ.r, [X.index...,UInt(m)])
    end
    return Y[1], Y[2:end]
end

function subdivide_unindexed(M::Integer, S::Vector{Similarity{T,M_}}, weights::Vector{Float64}, X::partition_data{T}) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
    Y = zeros(partition_data_unindexed{T},M)
    for (m,sₘ) in enumerate(S)
        Y[m] = partition_data_unindexed{T}(sim_map(sₘ,X.barycentre), X.weight*weights[m], X.diameter*sₘ.r)
    end
    return Y[1], Y[2:end]
end

subdivide(M::Integer,S::Vector{Similarity{T,M_}},weights::Vector{Float64}, X::partition_data_indexed{T}) where  {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}} = subdivide_indexed(M, S, weights, X)

subdivide(M::Integer, S::Vector{Similarity{T,M_}}, weights::Vector{Float64}, X::partition_data_unindexed{T}) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}} = subdivide_unindexed(M, S, weights, X)

function get_max_power(S::Vector{Similarity{T,M_}}, diameter::Float64, h::Float64) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
    r_max = 0.0
    for s in S r_max=max(r_max,s.r) end
    return ceil(Int64,log(h/diameter)/log(r_max))
end

function create_partition(γ::partition_data{T}, M::Int64, S::Vector{Similarity{T,M_}}, weights::Vector{Float64}, h::Float64) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
    max_size = M^get_max_power(S,γ.diameter, h)
    X = zeros(typeof(γ),max_size)
    X[1] = γ
    total_count = Int64(1)
    start_checking = Int64(1)
    subdiv = true
    while subdiv
        subdiv = false
        for n = start_checking:total_count
            if X[n].diameter > h
                range_extension = (total_count+1):(total_count+M-1)
                X[n], X[range_extension] = subdivide(M, S, weights, X[n])
                total_count += (M-1)
                subdiv = true
                break
            else
                start_checking += 1
                # break
            end
        end
    end
    return X[1:total_count]
end

function get_quadrature_from_partition(γ::partition_data{T}, M::Int64, S::Vector{Similarity{T,M_}}, weights::Vector{Float64}, h::Float64) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
    X = create_partition(γ, M, S, weights, h)
    N = length(X)
    w = zeros(N)
    zero_node = zero(γ.barycentre)
    x = [zero_node for _=1:N]
    for n=1:N
        x[n] = X[n].barycentre
        w[n] = X[n].weight
    end
    return x, w
end

barycentre_rule(Γ::Attractor{T,M_},h::Float64) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}} = get_quadrature_from_partition(partition_data_unindexed{T}(Γ.barycentre,Γ.measure,Γ.diameter), length(Γ.IFS), Γ.IFS, Γ.weights, h)

barycentre_rule(Γ::SubAttractor{T,M_},h::Float64) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}} = get_quadrature_from_partition(partition_data_unindexed{T}(Γ.barycentre,Γ.measure,Γ.diameter), length(Γ.Attractor.IFS), Γ.Attractor.IFS, Γ.Attractor.weights, h)

