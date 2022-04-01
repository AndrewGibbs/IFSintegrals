import Roots: find_zero, Bisection
import LinearAlgebra: norm, I, UniformScaling
import StaticArrays: SVector, SMatrix

abstract type ContractionMap
end

"""
Type Similarity(r,δ,A), with contraction r<1, translation δ ∈ Rⁿ and a rotation matrix A ∈ Rⁿˣⁿ.
"""
struct Similarity{V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    r::Float64 # contraction
    δ::V # translation
    A::M # rotation
    rA::M # contraction * rotation
end


function Similarity(r::Real, δ::Union{Vector{<:Real},Real}, θ=0::Real)
    # θ = Float64(θ)
    r = Float64(r)
    #build translation
    if isa(δ,Real)
        ndims = 1
        static_δ = Float64(δ)
    else
        ndims = length(δ)
        static_δ = SVector{ndims}(δ)
    end
    #build rotation
    if ndims ==1
        A = 1.0
    elseif θ == 0
        A = SMatrix{2,2,Float64}(I(2))
    else
        A = SMatrix{ndims,ndims}([cos(θ) -sin(θ); sin(θ) cos(θ)])
    end
    return Similarity(r,static_δ,A,r*A)
end

# # convert real values to vectors when appropriate
# Similarity(r::Real, δ::Real, θ=0::Real) = Similarity(r, [δ], θ)

#fast method, for N>2 and M=N^2 (although I can't work out how to restrict in this way)
# sim_map(s::Similarity{SVector{N,Float64},SMatrix{N,N,Float64,M}}, x::SVector{N,Float64}) where {N,M} = s.rA*x + s.δ

#fast method for N>2 without rotation
# sim_map(s::Similarity{SVector{N,Float64},UniformScaling{<:Real}}, x::SVector{N,Float64}) where N = s.r*x + s.δ

#fast method, for N=1
# sim_map(s::Similarity{Float64,Float64}, x::Float64) = s.r*x + s.δ

#abstract method, shouldn't really get used
sim_map(s::Similarity, x::Union{Real,AbstractVector{<:Real}}) = s.rA*x .+ s.δ

# sim_map(s::Similarity{1}, x::Float64) = s.r*x .+ s.δ[1]

fixed_point(s::Similarity{N}) where N = (I-s.r*s.A) \ s.δ
#s.δ/(1-s.r)

# there is possibility of defining affine contractions later.

# now define attractors of IFSs, as an subtype of a fractal
abstract type Fractal
end

# a couple of useful functions:
function get_dimension_from_contractions(R, homogeneous, top_dim)
    if homogeneous
        d = log(1/length(R))/log(R[1])
    else
        # approximate Hausdorff dimension:
        f(d) = sum(R.^d) - 1
        d = find_zero(f, (0,top_dim), Bisection())
    end
    return d
end

# function get_barycentre(sims::Array{Similarity{N}},d::Real) where N
#     M = length(sims[1].δ)
#     divisor = I
#     vec_sum = zeros(M)
#     for s in sims
#         divisor -= s.r^(d+1)*s.A
#         vec_sum += s.r^d*s.δ
#     end
#     return divisor \ vec_sum
# end

function get_barycentre(sims::Array{Similarity{N}}, weights::Vector{<:Real}) where N
    M = length(sims)
    divisor = I
    vec_sum = zeros(N)
    for n=1:M
        divisor -= sims[n].r*weights[n]*sims[n].A
        vec_sum += weights[n]*sims[n].δ
    end
    return divisor \ vec_sum
end


"""
    Attractor(sims::Array{Similarity}; measure::Real=1.0) = Attractor(sims, get_diameter(sims); measure=measure)
    
Representation of an attractor of an iterated function system (IFS).
Constructor requires only an IFS, which is of type Array{Similarity}, and diameter.
All other essential properties can be deduced from this, including barycentre
and dimension, which are approximated numerically.
"""
struct Attractor{M,N} <: Fractal
    IFS::Vector{Similarity{N}}
    topological_dimension::Int64
    Hausdorff_dimension::Float64
    Hausdorff_weights::Bool
    homogeneous::Bool
    barycentre::SVector{N,Float64}
    diameter::Float64
    measure::Float64
    weights::SVector{M,Float64}
end

# outer constructor, when diameter isn't given:
"""
This rests on the assumption that the convex hull of the fractal
is the convex hull of its fixed points. I haven't been able to 
construct an example where this is false yet, but I haven't considered
any rotations.
"""

# Attractor(sims::Array{Similarity}; measure::Real=1.0) = Attractor(sims, get_diameter(sims); measure=measure)

function Attractor(sims::Vector{Similarity{N}}; diameter::Real=0.0, measure::Real=1.0, weights::Vector{<:Real}=[0]) where N# outer constructor for attractor type
    count = 1
    top_dims = zeros(Int64,length(sims))
    contractions = zeros(Float64,length(sims))
    for s in sims
        if s.r>=1
            error("Contraction factor must be <1")
        end
        top_dims[count] = length(s.δ)
        contractions[count] = s.r
        count += 1
    end
    if !constant_vector(top_dims)
        error("Contractions are of different topological dimension")
    end

    if weights != [0]
        if length(weights)!=length(sims)
            error("Weights vector must be the same length as similarities vector")
        end
    end

    # Now construct the properties needed for manipulations which will follow.

    # Topological dimension
    top_dim = maximum(top_dims)
    if constant_vector(contractions)
        homogeneous = true
    else
        homogeneous = false
    end

    #Hausdorff dimension
    Hdim = get_dimension_from_contractions(contractions,homogeneous,top_dim)

    #Barycentre
    bary = get_barycentre(sims,weights)
    Sbary = SVector{top_dim,Float64}(bary)

    #Diameter
    if diameter <= 0.0
        diameter = get_diameter(sims)
    end

    M = length(sims)
    if weights == [0]
        weights = zeros(Float64,M)
        for m=1:M
            weights[m] = measure*sims[m].r^Hdim
        end
    end
    Sweights = SVector{M,Float64}(weights)

    Hausdorff_weights = are_weights_Hausdorff(weights, sims, Hdim)

    return Attractor{M,top_dim}(sims, top_dim, Hdim, homogeneous, Hausdorff_weights, Sbary, Float64(diameter), Float64(measure), Sweights)

end

# subcomponent of attractor, as a subclass of fractal
struct SubAttractor{M,N} <: Fractal
    attractor::Attractor{M,N}
    # IFS::Array{Similarity} # could be removed
    index::Vector{UInt8}
    # topological_dimension::Int64 # could be removed
    # Hausdorff_dimension::T # could be removed
    barycentre::SVector{N,Float64}
    diameter::Float64
    measure::Float64
end
# outer constructor
"""
Representation of a subcomponent of a fractal Γ, using standard vector index notatation.
If Γ is a subattractor, then the vector indices are concatenated to produce a new subatractor,
which stores the original attractor.
"""
function SubAttractor(Γ::Union{Attractor{M,N},SubAttractor{M,N}}, index::Vector{UInt8}) where {M,N}

    #quick condition for trivial case:
    if index == [0]
        if typeof(Γ)==SubAttractor
            return Γ
        else
            return SubAttractor(Γ, [UInt8(0)], Γ.barycentre, Γ.diameter, Γ.measure, Γ.homogeneous)
        end
    else #non-trivial case

        # get new measure and diameter. First initialise:
        new_diam = Γ.diameter
        new_measure = Γ.measure
        # now adjust the diameter and measure for the smaller scale:
        for m = index
            new_diam *= Γ.IFS[m].r
            new_measure *= Γ.weights[m]
        end

        # get new vector index and barycentre:
        if typeof(Γ)==SubAttractor
            if Γ.index != [0] # excluding trivial case
                index = vcat(Γ.index,index)
            end
            orig_bary = Γ.attractor.barycentre
        else
            orig_bary = Γ.barycentre
        end

        #start as old barycentre and map
        new_bary = orig_bary
        for m=index[end:-1:1]
            new_bary = sim_map(Γ.IFS[m], new_bary)
        end
        Snew_bary = SVector{N,Float64}(new_bary)

        #quick fix if we're after a subattractor of a subattractor
        if typeof(Γ)==SubAttractor
            #index = vcat(index, Γ.index)
            return SubAttractor{M,N}(Γ.attractor, index, Snew_bary, new_diam, new_measure)
        else
            return SubAttractor{M,N}(Γ,           index, Snew_bary, new_diam, new_measure)
        end
    end
end

SubAttractor(Γ::Union{Attractor,SubAttractor}, index::Int64) = SubAttractor(Γ, [index])

full_map(S::Array{Similarity{N}},x::Array{<:Real,1}) where N = full_map(S,reshape(x,length(x),1))

function are_weights_Hausdorff(w::Vector{Float64},S::Vector{Similarity{N}},d::Number) where N
    thresh = 1E-12
    x = true
    for n=1:length(w)
        if abs(w[n]-S[n].r^d) > thresh
            x = false
            break
        end
    end
    return x
end

function full_map(S::Vector{Similarity{N}},X::Vector{Vector{Float64}}) where N
    d = length(X)
    M = length(S)
    Y = [Vector{Float64}(undef, N) for _ = 1:(M*d)]
    for nₓ = 1:d
        for m=1:M
            Y[d*(m-1)+nₓ] = sim_map(S[m], X[nₓ])
        end
    end
    return Y
end

# provides a sketch of an attractor in N topological dimensions
function sketch_attractor(γ::Attractor; mem_const = 10000)
    N = 10
    X = rand(γ.topological_dimension,N)
    M = length(γ.IFS)
    num_its = floor(log(mem_const/(γ.topological_dimension*N))/log(M))
    #println(num_its)
    S = γ.IFS
    count = 0
    for count = 1:num_its
        NxM = N*M
        X_ = zeros(γ.topological_dimension,NxM)
        s_count = 0
        for m=1:M
            X_[:,((m-1)*N+1):(m*N)] = sim_map(S[m], X)
        end
        X = X_
        N = NxM
    end
    return X
end

# now define attractors of popular fractals
"""
The middle-α Cantor set (default is α=1/3),
formed by removing the middle α of the unit interval, and repeating on each sub interval.
"""

# quick function which computes the barycentre, if it's not the default one with equal weights
function default_bary(S::Vector{Similarity{N}}, d::Number, weights::Vector{Float64}, Hausdorff_measure_bary::Vector{Float64}) where N
    if are_weights_Hausdorff(weights,S,d)
        bary = Hausdorff_measure_bary
    else
        bary = get_barycentre(S,weights)
    end
    return bary
end

function CantorSet(;contraction = 1/3, weights=[1/2, 1/2])
    S = [Similarity(contraction,[0.0]),Similarity(contraction,[1-contraction])]
    d = log(1/2)/log(contraction)
    bary = default_bary(S,d,weights,[1/2])
    return Attractor{2,1}(S,1,d,true,are_weights_Hausdorff(weights,S,d),bary,1.0,1.0,weights)
end

"""
The middle-α Cantor dust (default is α=1/3),
formed by taking the cartesian product of two middle-α Cantor sets.
"""
function CantorDust(;contraction = 1/3, weights=[1/4, 1/4, 1/4, 1/4])
    S = [Similarity(contraction,[0.0,0.0]),Similarity(contraction,[1-contraction,0.0]),Similarity(contraction,[0.0,1-contraction]),Similarity(contraction,[1-contraction,1-contraction])]
    d = log(1/4)/log(contraction)
    bary = default_bary(S,d,weights,[0.5,0.5])
    return Attractor{4,2}(S, 2, d, true, are_weights_Hausdorff(weights,S,d), bary, sqrt(2), 1.0,weights)
end

function CantorN(N::Integer; contraction = 1/3)
    M = 2^N
    S = Similarity[]
    for m=1:M
        m_binary = digits(m-1,base=2,pad=N)
        push!(S,Similarity(contraction,m_binary.*(1-contraction)))
    end
    d = log(1/M)/log(contraction)
    weights = ones(M)./M
    bary =  0.5.*ones(N)
    return Attractor{M,N}(S, 2, d, true, true, bary, sqrt(N), 1.0, weights)
end

"""
The Sierpinski triangle, as an attractor of an iterated function system.
"""
function Sierpinski(;weights=[1/3, 1/3, 1/3])
    courage = Similarity(1/2,[0,1/6])
    wisdom = Similarity(1/2,sqrt(2)*[1/6,-1/6])
    power = Similarity(1/2,sqrt(2)*[-1/6,-1/6])
    S = [courage,wisdom,power]
    d = log(3)/log(2)
    bary = default_bary(S,d,weights,[0,(1-2*sqrt(2))/9])
    return Attractor{3,2}(S, 2, d, true, are_weights_Hausdorff(weights,S,d), bary, 1.0, 1.0, weights)
end
#            Attractor(sims,top_dim,Hdim,uniform,get_barycentre(sims,Hdim),diameter,measure)

function SquareFlake(;weights=ones(16)./16)
    scale = 1
    h = scale/4 # width of square subcomponent
    ρ = 1/4
    IFS = [ Similarity(ρ,[-3h/2,3h/2]), #1
            Similarity(ρ,[-h/2,5h/2]), #2
            Similarity(ρ,[-h/2,3h/2]), #3
            Similarity(ρ,[3h/2,3h/2]), #4
            Similarity(ρ,[-h/2,h/2]), #5
            Similarity(ρ,[h/2,h/2]), #6
            Similarity(ρ,[3h/2,h/2]), #7
            Similarity(ρ,[5h/2,h/2]), #8
            Similarity(ρ,[-5h/2,-h/2]), #9
            Similarity(ρ,[-3h/2,-h/2]), #10
            Similarity(ρ,[-h/2,-h/2]), #11
            Similarity(ρ,[h/2,-h/2]), #12
            Similarity(ρ,[-3h/2,-3h/2]), #13
            Similarity(ρ,[h/2,-3h/2]), #14
            Similarity(ρ,[3h/2,-3h/2]), #15
            Similarity(ρ,[h/2,-5h/2]) #16
            ]
    R = get_diameter(IFS) # I'm sure this can be calculated by hand... but not today.
    # also the 'measure' is not really 1 here. But it doesn't matter.
    bary = default_bary(IFS,2.0,weights,[0.0,0.0])
    return Attractor{16,2}(IFS, 2, 2.0, true, are_weights_Hausdorff(weights,IFS,2), [0.0,0.0], R, 1.0, weights)
end

function KochFlake(weights = [1/3, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9])
    IFS = [Similarity(sqrt(1/3),[0, 0], pi/6),
            Similarity(1/3,[1/sqrt(3),1/3]),
            Similarity(1/3,[0,2/3]),
            Similarity(1/3,[-1/sqrt(3),1/3]),
            Similarity(1/3,[-1/sqrt(3),-1/3]),
            Similarity(1/3,[0,-2/3]),
            Similarity(1/3,[1/sqrt(3),-1/3])
            ]
    bary = default_bary(IFS,2.0,weights,[0.0,0.0])
    return Attractor{7,2}(IFS, 2, 2.0, false, are_weights_Hausdorff(weights,IFS,2), bary, 2*sqrt(3)/3, 1.0, weights)
end