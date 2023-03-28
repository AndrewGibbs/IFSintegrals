"""
    struct Similarity{V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
        r::Float64 # contraction
        δ::V # translation
        A::M # rotation/reflection
        rA::M # contraction * rotation/reflection
    end

Constructs a similarity map. 
The third input (rotation matrix) is optional, and the fourth is created automatically.


Can be treated as a function, for example
```julia-repl
    julia> s = Similarity(1/2,0) # creates contraction s(x)=x/2+0
    julia> s(π)
    1.5707963267948966
```
"""
struct Similarity{V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    r::Float64 # contraction
    δ::V # translation
    A::M # rotation
    rA::M # contraction * rotation
end

# constructor
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

# define similarity as a map
sim_map(s::Similarity, x) = s.rA*x + s.δ
(s::Similarity{V,M})(x::V) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}} = sim_map(s,x)
#inverse map
sim_map_inv(s::Similarity, x) = s.rA/(x-s.δ)

# similarity acting on an IFS
function sim_map(s::Similarity{V,M}, IFS::Vector{Similarity{V,M}}) where {V<:Union{Real,AbstractVector}, M<:AbstractMatrix}
    Achain = [s.A*IFS[m].A/s.A for m=1:length(IFS)]
    rAchain = [IFS[m].r*Achain[m] for m=1:length(IFS)]
    return [Similarity{V,M}(IFS[m].r, (I-rAchain[m])*s.δ + s.rA*IFS[m].δ, Achain[m], rAchain[m]) for m=1:length(IFS)]
end
# simplified case without rotation
function sim_map(s::Similarity{V,M}, IFS::Vector{Similarity{V,M}}) where {V<:Union{Real,AbstractVector}, M<:Real}
    return [Similarity{V,M}(IFS[m].r, (I-IFS[m].rA)*s.δ + s.rA*IFS[m].δ, IFS[m].A, IFS[m].rA) for m=1:length(IFS)]
end

# composition of two similarities
∘(s₁::Similarity{V,M},s₂::Similarity{V,M}) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}} = Similarity(s₁.r*s₂.r, s₁.δ+s₁.rA*s₂.δ, s₁.A*s₂.A, s₁.rA*s₂.rA)

function sim_comp(IFS::Vector{Similarity{V,M}}, m::Vector{Int64}) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    s = IFS[m[1]]
        for j = 2:length(m)
            s = s∘IFS[m[j]]
        end
    return s
end

# union of similarities acting on a set of points, returning a set of points:
function full_map(S::Vector{Similarity{V,M_}},X::Vector{V_}) where {V<:AbstractVector{<:Real}, V_<:AbstractVector{<:Real}, M_<:Union{Real,AbstractMatrix}}
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
function full_map(S::Vector{Similarity{V,M_}},X::Vector{V_}) where {V<:Real, V_<:Real, M_<:Real}
    d = length(X)
    M = length(S)
    Y = zeros(V_,d*M)
    for nₓ = 1:d
        for m=1:M
            Y[d*(m-1)+nₓ] = sim_map(S[m], X[nₓ])
        end
    end
    return Y
end

# versions acting on a single vector/value
full_map(S::Vector{Similarity{V,M}},X::V_) where {V<:AbstractVector, V_<:AbstractVector, M<:Union{Real,AbstractMatrix}} = (full_map(S,[X]))
full_map(S::Vector{Similarity{V,M}},X::V_) where {V<:Real, V_<:Real, M<:Real} = (full_map(S,[X]))

# shorthand for everything, which appears quite natural
(S::Vector{Similarity{V,M_}})(X) where {V<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}} = full_map(S,X)

# define the fixed point of the map
fixed_point(s::Similarity{V,M}) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}} = (I-s.rA) \ s.δ