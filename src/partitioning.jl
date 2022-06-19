import Base: zero

abstract type partition_data{T<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}}
end

struct partition_data_indexed{T<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}} <: partition_data{T,M}
    barycentre::T
    weight::Float64
    diameter::Float64
    index::Vector{Int64}
end

struct partition_data_unindexed{T<:Union{Real,AbstractVector},M<:Union{Real,AbstractMatrix}} <: partition_data{T,M}
    barycentre::T
    weight::Float64
    diameter::Float64
end

struct partition_data_with_IFS{T<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}} <: partition_data{T,M}
    barycentre::T
    weight::Float64
    diameter::Float64
    IFS::Vector{Similarity{T,M}}
end

zero(::Type{partition_data_indexed{T}}) where {T<:Union{Real,AbstractVector}} =  partition_data_indexed{T}(zero(T),0.0,0.0,[0])

zero(::Type{partition_data_unindexed{T}}) where {T<:Union{Real,AbstractVector}} =  partition_data_unindexed{T}(zero(T),0.0,0.0)

function subdivide(M::Integer,S::Vector{Similarity{T,M_}},weights::Vector{Float64}, X::partition_data_indexed{T}) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
    Y = zeros(partition_data_indexed{T},M)
    for m=1:M
        Y[m] = partition_data_indexed{T}(sim_map(S[m],X.barycentre), X.weight*weights[m], X.diameter*S[m].r, [X.index...,m])
    end
    return Y[1], Y[2:end]
end

function subdivide(M::Integer, S::Vector{Similarity{T,M_}}, weights::Vector{Float64}, X::partition_data_unindexed{T}) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
    Y = zeros(partition_data_unindexed{T},M)
    for m=1:M
        Y[m] = partition_data_unindexed{T}(sim_map(S[m],X.barycentre), X.weight*weights[m], X.diameter*S[m].r)
    end
    return Y[1], Y[2:end]
end

function subdivide(M::Integer, S::Vector{Similarity{T,M_}}, weights::Vector{Float64}, X::partition_data_with_IFS{T,M_}) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
    IFSnew = [sim_map(S[m],X.IFS) for m=1:M]
    Y =[partition_data_with_IFS{T,M_}(sim_map(S[m],X.barycentre), X.weight*weights[m], X.diameter*S[m].r, IFSnew[m]) for m=1:M]
    return Y[1], Y[2:end]
end

# subdivide(M::Integer,S::Vector{Similarity{T,M_}},weights::Vector{Float64}, X::partition_data_indexed{T}) where  {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}} = subdivide_indexed(M, S, weights, X)

# subdivide(M::Integer, S::Vector{Similarity{T,M_}}, weights::Vector{Float64}, X::partition_data_unindexed{T}) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}} = subdivide_unindexed(M, S, weights, X)

function get_max_power(S::Vector{Similarity{T,M_}}, diameter::Float64, h::Float64) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
    r_max = 0.0
    for s in S r_max=max(r_max,s.r) end
    return max(ceil(Int64,log(h/diameter)/log(r_max)),0)
end

function create_partition(γ::partition_data{T}, M::Int64, S::Vector{Similarity{T,M_}}, weights::Vector{Float64}, h::Float64) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
    max_size = M^get_max_power(S,γ.diameter*(1+eps()), h)
    # X = zeros(typeof(γ),max_size)
    X = [γ for _=1:max_size]
    # X[1] = γ
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

function subdivide_indices(Γ::SelfSimilarFractal, h::Real; int_type::DataType=Int64)
    I = Vector{int_type}[]
    if isa(Γ,Attractor)
        IFS = Γ.IFS
    else
        IFS = Γ.attractor.IFS
    end
    M = length(IFS)
    r = zeros(M)

    if Γ.diameter >= h
        subdiv = true
        for m=int_type.(1:M)
            push!(I,[m])
            r[m] = IFS[m].r
        end
    else
        subdiv = false
    end

    while subdiv
        subdiv = false
        split_vecs = Int64[]
        for j = 1:length(I)
           if Γ.diameter*prod(r[I[j]]) >= h
                subdiv = true
                push!(split_vecs,j)
            end
        end
        if subdiv
            new_vecs = Vector{int_type}[]
            for j in split_vecs
                for m = 1:M
                    push!(new_vecs,vcat(I[j],[m]))
                end
            end
            deleteat!(I,split_vecs)
            I = vcat(I,new_vecs)
        end
    end
    #quick bodge - this convention means we can keep the same type
        # and it's (more) consistent with the paper
    if isempty(I)
        I = [[int_type(0)]]
    end
    return I#convert(Array{Array{Int64,1},1},I)
end