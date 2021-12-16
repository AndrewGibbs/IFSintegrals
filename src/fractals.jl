import Roots: find_zero, Bisection
import LinearAlgebra: norm, I, UniformScaling

abstract type ContractionMap
end

"""
Type Similarity(r,δ,A), with contraction r<1, translation δ ∈ Rⁿ and a rotation matrix A ∈ Rⁿˣⁿ.
"""
struct Similarity <: ContractionMap
    r::Real # contraction
    δ::Array{<:Real,1} # translation
    A::Union{Array{<:Real,2},UniformScaling{Bool}} # rotation
    Similarity(r,δ,A) = r >= 1 || r<=0 ? error("Contraction factor must be ∈(0,1)") : new(r,δ,A)
    Similarity(r,δ,A) = typeof(A)!=UniformScaling{Bool} ? (norm(A*A'-I) > 1E-8 ? error("Rotation matrix must be orthoganol") : new(r,δ,A)) : new(r,δ,A)
end

# Similarity(r::Real, δ::Union{Real,Array{<:Real,1}}) = Similarity(r,δ,Matrix(1.0I,length(δ),length(δ)))

# Similarity(r::Real, δ::Real, A::Array{<:Real,2}) = Similarity(r,[δ],A)

function Similarity(r::Real, δ::Array{<:Real,1}; θ=0::Real)
    if θ == 0
        return Similarity(r,δ,I)
    else
        A = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        return Similarity(r,δ,A)
    end
end

sim_map(s::Similarity, x::Union{Array{<:Real},AbstractArray{<:Real,2}}) = s.r*s.A*x .+ s.δ

fixed_point(s::Similarity) = (I-s.r*s.A) \ s.δ
#s.δ/(1-s.r)

# there is possibility of defining affine contractions later.

# now define attractors of IFSs, as an subtype of a fractal
abstract type Fractal
end

# a couple of useful functions:
function get_dimension_from_contractions(R, uniform, top_dim)
    if uniform
        d = log(1/length(R))/log(R[1])
    else
        # approximate Hausdorff dimension:
        f(d) = sum(R.^d) - 1
        d = find_zero(f, (0,top_dim), Bisection())
    end
    return d
end

function get_barycentre(sims::Array{Similarity},d::Real)
    N = length(sims[1].δ)
    divisor = I
    vec_sum = zeros(N)
    for s in sims
        divisor -= s.r^(d+1)*s.A
        vec_sum += s.r^d*s.δ
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
struct Attractor <: Fractal
    IFS::Array{Similarity}
    topological_dimension::Int64
    Hausdorff_dimension::Real
    uniform::Bool
    barycentre::Array{<:Real}
    diameter::Real
    measure::Real
end

# outer constructor, when diameter isn't given:
"""
This rests on the assumption that the convex hull of the fractal
is the convex hull of its fixed points. I haven't been able to 
construct an example where this is false yet, but I haven't considered
any rotations.
"""

Attractor(sims::Array{Similarity}; measure::Real=1.0) = Attractor(sims, get_diameter(sims); measure=measure)

function Attractor(sims::Array{Similarity}, diameter::Real; measure::Real=1.0) # outer constructor for attractor type
    count = 1
    top_dims = zeros(length(sims))
    contractions = zeros(length(sims))
    for s in sims
        if s.r>=1
            error("Contraction factor must be <1")
        end
        top_dims[count] = length(s.δ)
        contractions[count] = s.r
        count += 1
    end
    if maximum(top_dims) != minimum(top_dims)
        error("Contractions are of different topological dimension")
    else
        top_dim = maximum(top_dims)
        if maximum(contractions) != minimum(contractions)
            uniform = false
        else
            uniform = true
        end
        Hdim = get_dimension_from_contractions(contractions,uniform,top_dim)
        return Attractor(sims,top_dim,Hdim,uniform,get_barycentre(sims,Hdim),diameter,measure)
    end
end

# subcomponent of attractor, as a subclass of fractal
struct SubAttractor <: Fractal
    attractor::Attractor
    IFS::Array{Similarity}
    index::Array{Int64}
    topological_dimension::Int64 # could be removed
    Hausdorff_dimension::Real # could be removed
    barycentre::Array{<:Real}
    diameter::Real
    measure::Real
    uniform::Bool
end
# outer constructor
"""
Representation of a subcomponent of a fractal Γ, using standard vector index notatation.
If Γ is a subattractor, then the vector indices are concatenated to produce a new subatractor,
which stores the original attractor.
"""
function SubAttractor(Γ::Union{Attractor,SubAttractor}, index::Array{Int64})

    #quick condition for trivial case:
    if index == [0]
        if typeof(Γ)==SubAttractor
            return Γ
        else
            return SubAttractor(Γ, Γ.IFS, [0], Γ.topological_dimension,
                            Γ.Hausdorff_dimension, Γ.barycentre,
                            Γ.diameter, Γ.measure, Γ.uniform)
        end
    else #non-trivial case

        # get new measure and diameter:
        R = 1
        for m=index
            R *= Γ.IFS[m].r
        end
        new_diam = R * Γ.diameter
        new_measure = R^(Γ.Hausdorff_dimension) * Γ.measure

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

        #quick fix if we're after a subattractor of a subattractor
        if typeof(Γ)==SubAttractor
            #index = vcat(index, Γ.index)
            return SubAttractor(Γ.attractor, Γ.IFS, index, Γ.topological_dimension, Γ.Hausdorff_dimension, new_bary, new_diam, new_measure, Γ.uniform)
        else
            return SubAttractor(Γ,           Γ.IFS, index, Γ.topological_dimension, Γ.Hausdorff_dimension, new_bary, new_diam, new_measure, Γ.uniform)
        end
    end
end
SubAttractor(Γ::Union{Attractor,SubAttractor}, index::Int64) = SubAttractor(Γ, [index])

full_map(S::Array{Similarity},x::Array{<:Real,1}) = full_map(S,reshape(x,length(x),1))

function full_map(S::Array{Similarity},X::Array{<:Real,2})
    d,N = size(X)
    M = length(S)
    NxM = N*M
    Y = zeros(d,NxM)
    for m=1:M
        Y[:,((m-1)*N+1):(m*N)] = sim_map(S[m], X)
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
function CantorSet(α = 1/3)
    r = (1 - α)/2
    S = [Similarity(r,[0.0]),Similarity(r,[1-r])]
    return Attractor(S,1,log(1/2)/log(r),true,[1/2],1.0,1.0)
end

"""
The middle-α Cantor dust (default is α=1/3),
formed by taking the cartesian product of two middle-α Cantor sets.
"""
function CantorDust(α = 1/3)
    r = (1 - α)/2
    S = [Similarity(r,[0.0,0.0]),Similarity(r,[1-r,0.0]),Similarity(r,[0.0,1-r]),Similarity(r,[1-r,1-r])]
    return Attractor(S, 2, log(1/4)/log(r), true, [0.5,0.5], sqrt(2), 1.0)
end

"""
The Sierpinski triangle, as an attractor of an iterated function system.
"""
function Sierpinski()
    courage = Similarity(1/2,[0,1/6])
    wisdom = Similarity(1/2,sqrt(2)*[1/6,-1/6])
    power = Similarity(1/2,sqrt(2)*[-1/6,-1/6])
    return Attractor([courage,wisdom,power], 2, log(3)/log(2), true, [0,(1-2*sqrt(2))/9], 1.0, 1.0)
end
#            Attractor(sims,top_dim,Hdim,uniform,get_barycentre(sims,Hdim),diameter,measure)

function SquareFlake()
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
    return Attractor(IFS, 2, 2, true, [0.0,0.0], R, 1.0)
end

function KochFlake()
    IFS = [Similarity(sqrt(1/3),[0, 0];θ = pi/6),
            Similarity(1/3,[1/sqrt(3),1/3]),
            Similarity(1/3,[0,2/3]),
            Similarity(1/3,[-1/sqrt(3),1/3]),
            Similarity(1/3,[-1/sqrt(3),-1/3]),
            Similarity(1/3,[0,-2/3]),
            Similarity(1/3,[1/sqrt(3),-1/3])
            ]
    return Attractor(IFS, 2, 2, false, [0.0,0.0], 2*sqrt(3)/3, 1.0)
end