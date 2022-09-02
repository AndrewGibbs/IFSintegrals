struct Similarity{V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    r::Float64 # contraction
    δ::V # translation
    A::M # rotation
    rA::M # contraction * rotation
end

function Similarity(r::Real, δ::Union{Vector{<:Real},Real}, rotation::Union{AbstractMatrix,Real}=I(length(δ)))
    # θ = Float64(θ)
    r = Float64(r)
    #build translation
    if isa(δ,Real)
        ndims = 1
        static_δ = Float64(δ)
    else
        ndims = length(δ)
        static_δ = SVector{ndims,Float64}(δ)
    end

    if isa(rotation,Real) # convert to rotation matrix
        if ndims != 2
            error("If rotation is to be determined from a single angle, similarity must be two-dimensional")
        end
        rotation = [cos(rotation) -sin(rotation); sin(rotation) cos(rotation)]
    end

    #build rotation
    if ndims ==1
        A = 1.0
    else
        A = SMatrix{ndims,ndims,Float64}(rotation)
    end
    return Similarity(r,static_δ,A,r*A)
end

# function Similarity(r::Real, δ::Union{Vector{<:Real},Real}, θ=0::Real)
#     ndims = length(δ)
#     if θ == 0
#         return Similarity(r, δ, I(ndims))
#     else
#         if ndims != 2
#             error("If rotation is to be determined from a single angle, similarity must be two-dimensional")
#         end
#         A = [cos(θ) -sin(θ); sin(θ) cos(θ)]
#         return Similarity(r, δ, A)
#     end
# end

sim_map(s::Similarity, x::Union{Real,AbstractVector{<:Real}}) = s.rA*x + s.δ
sim_map_inv(s::Similarity, x::Union{Real,AbstractVector{<:Real}}) = s.rA/(x-s.δ)

function sim_map(s::Similarity{V,M}, IFS::Vector{Similarity{V,M}}) where {V<:Union{Real,AbstractVector}, M<:AbstractMatrix}
    Achain = [s.A*IFS[m].A/s.A for m=1:length(IFS)]
    rAchain = [IFS[m].r*Achain[m] for m=1:length(IFS)]
    return [Similarity{V,M}(IFS[m].r, (I-rAchain[m])*s.δ + s.rA*IFS[m].δ, Achain[m], rAchain[m]) for m=1:length(IFS)]
end

∘(s₁::Similarity{V,M},s₂::Similarity{V,M}) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}} = Similarity(s₁.r*s₂.r, s₁.δ+s₁.rA*s₂.δ, s₁.A*s₂.A, s₁.rA*s₂.rA)
"""   
r::Float64 # contraction
δ::V # translation
A::M # rotation
rA::M # contraction * rotation
"""
function sim_comp(IFS::Vector{Similarity{V,M}}, m::Vector{Int64}) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    s = IFS[m[1]]
        for j = 2:length(m)
            s = s∘IFS[m[j]]
        end
    return s
end

# below could just be rs
function sim_map(s::Similarity{V,M}, IFS::Vector{Similarity{V,M}}) where {V<:Union{Real,AbstractVector}, M<:Real}
    return [Similarity{V,M}(IFS[m].r, (I-IFS[m].rA)*s.δ + s.rA*IFS[m].δ, IFS[m].A, IFS[m].rA) for m=1:length(IFS)]
end

function full_map(S::Vector{Similarity{V,M_}},X::Vector{V_}) where {V<:Union{Real,AbstractVector}, V_<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
    d = length(X)
    M = length(S)
    N = length(S[1].δ)
    Y = [0*X[1] for _ = 1:(M*d)]
    for nₓ = 1:d
        for m=1:M
            Y[d*(m-1)+nₓ] = sim_map(S[m], X[nₓ])
        end
    end
    return Y
end

full_map(S::Vector{Similarity{V,M_}},X::V) where {V<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}} = (full_map(S,[X]))[1]

# full_map(S::Array{Similarity},x::Vector{<:Real}) where N = full_map(S,reshape(x,length(x),1))

fixed_point(s::Similarity{V,M}) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}} = (I-s.rA) \ s.δ