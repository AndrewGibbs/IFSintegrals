"""
Type Similarity(r,δ,A), with contraction r<1, translation δ ∈ Rⁿ and a rotation matrix A ∈ Rⁿˣⁿ.
"""

# a few useful functions for fractal properties:
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

function get_barycentre(sims::Array{Similarity{V,M_}}, weights::Vector{<:Real}) where {V<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
    M = length(sims)
    N = length(sims[1].δ)
    divisor = I
    vec_sum = zeros(N)
    for n=1:M
        divisor -= sims[n].r*weights[n]*sims[n].A
        vec_sum += weights[n]*sims[n].δ
    end
    return divisor \ vec_sum
end

function are_weights_Hausdorff(w::Vector{Float64},S::Vector{Similarity{V,M}},d::Number) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
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


abstract type SelfSimilarFractal{V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
end

"""
    Attractor(sims::Array{Similarity}; measure::Real=1.0) = Attractor(sims, get_diameter(sims); measure=measure)
    
Representation of an attractor of an iterated function system (IFS).
Constructor requires only an IFS, which is of type Array{Similarity}, and diameter.
All other essential properties can be deduced from this, including barycentre
and dimension, which are approximated numerically.
"""
struct Attractor{V,M} <: SelfSimilarFractal{V,M}
    IFS::Vector{Similarity{V,M}}
    topological_dimension::Int64
    Hausdorff_dimension::Float64
    Hausdorff_weights::Bool
    homogeneous::Bool
    barycentre::V
    diameter::Float64
    measure::Float64
    weights::Vector{Float64}
end

# outer constructor, when diameter isn't given:
"""
This rests on the assumption that the convex hull of the fractal
is the convex hull of its fixed points. I haven't been able to 
construct an example where this is false yet, but I haven't considered
any rotations.
"""

# Attractor(sims::Array{Similarity}; measure::Real=1.0) = Attractor(sims, get_diameter(sims); measure=measure)

function Attractor(sims::Vector{Similarity{V,M_}}; diameter::Real=0.0, measure::Real=1.0, weights::Vector{<:Real}=[0.0]) where {V<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}# outer constructor for attractor type
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

    # weights 
    M = length(sims)
    if weights == [0]
        weights = zeros(Float64,M)
        for m=1:M
            weights[m] = measure*sims[m].r^Hdim
        end
    end
    # Sweights = SVector{M,Float64}(weights)

    Hausdorff_weights = are_weights_Hausdorff(weights, sims, Hdim)

    #Barycentre
    Sbary = get_barycentre(sims,weights)
    # Sbary = SVector{top_dim,Float64}(bary)

    #Diameter
    if diameter <= 0.0
        diameter = get_diameter(sims)
    end

    return Attractor(sims, top_dim, Hdim, homogeneous, Hausdorff_weights, Sbary, Float64(diameter), Float64(measure), weights)

end

# subcomponent of attractor, as a subclass of fractal
struct SubAttractor{V,M} <: SelfSimilarFractal{V,M}
    attractor::Attractor
    # IFS::Array{Similarity} # could be removed ?
    index::Vector{UInt64}
    # topological_dimension::Int64 # could be removed
    # Hausdorff_dimension::T # could be removed
    barycentre::Vector{Float64}
    diameter::Float64
    measure::Float64
end
# outer constructor
"""
Representation of a subcomponent of a fractal Γ, using standard vector index notatation.
If Γ is a subattractor, then the vector indices are concatenated to produce a new subatractor,
which stores the original attractor.
"""
function SubAttractor(Γ::Attractor{V,M}, index::Vector{<:Unsigned}) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    #quick condition for trivial case:
    if index == [0]
        return SubAttractor(Γ, [0], Γ.barycentre, Γ.diameter, Γ.measure)
    else #non-trivial case

        # get new measure and diameter. First initialise:
        new_diam = Γ.diameter
        new_measure = Γ.measure
        # now adjust the diameter and measure for the smaller scale:
        for m = index
            new_diam *= Γ.IFS[m].r
            new_measure *= Γ.weights[m]
        end

        #start as old barycentre and map
        new_bary = Γ.barycentre
        for m=index[end:-1:1]
            new_bary = sim_map(Γ.IFS[m], new_bary)
        end

        return SubAttractor{V,M}(Γ,           index, new_bary, new_diam, new_measure)
        
    end
end

function SubAttractor(Γ::SubAttractor{V,M}, index::Vector{<:Unsigned}) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
        #quick condition for trivial case:
        if index == [0]
            return Γ
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
            if Γ.index != [0] # excluding trivial case
                index = vcat(Γ.index,index)
            end
    
            #start as old barycentre and map
            new_bary = Γ.attractor.barycentre
            for m=index[end:-1:1]
                new_bary = sim_map(Γ.IFS[m], new_bary)
            end
    
            return SubAttractor{V,M}(Γ.attractor, index, new_bary, new_diam, new_measure)
        end
    end

SubAttractor(Γ::SelfSimilarFractal, index::Unsigned) = SubAttractor(Γ, [index])

# now define attractors of popular fractals
"""
The middle-α Cantor set (default is α=1/3),
formed by removing the middle α of the unit interval, and repeating on each sub interval.
"""

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