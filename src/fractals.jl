"""
Type Similarity(r,δ,A), with contraction r<1, translation δ ∈ Rⁿ and a rotation matrix A ∈ Rⁿˣⁿ.
"""

# a few useful functions for fractal properties:
function get_dimension_from_contractions(R, homogeneous, top_dim)
    if sum(R.^top_dim)>1+100*eps()
        error("Contractions are too large - resulting fractal will have Hausdorff dimension larger than spatial dimension.")
    end
    if homogeneous
        d = log(1/length(R))/log(R[1])
    else
        # approximate Hausdorff dimension:
        f(d) = sum(R.^d) - 1
        d = find_zero(f, (0,top_dim*(1+eps())), Bisection())
        # have added +ϵ here, to ensure that endpoints of search are of opposite sign when d=n
    end
    return d
end

function get_barycentre(sims::Array{Similarity{V,M_}}, weights::Vector{<:Real}) where {V<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
    M = length(sims)
    divisor = I
    if V<:Real
        vec_sum = 0.0
    else
        vec_sum = zero(sims[1].δ)
    end
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


abstract type FractalMeasure{V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}} end

abstract type SelfSimilarFractal{V,M}<:FractalMeasure{V,M} end

"""
    struct InvariantMeasure{V,M} <: SelfSimilarFractal{V,M}
        IFS::Vector{Similarity{V,M}}
        spatial_dimension::Int64
        Hausdorff_dimension::Float64
        homogeneous::Bool
        Hausdorff_weights::Bool
        barycentre::V
        diameter::Float64
        measure::Float64
        weights::Vector{Float64}
        disjoint::Bool
        connectedness::Matrix{Bool}
        symmetry_group::Vector{AutomorphicMap{V,M}}
    end
    
Representation of a invariant measure, whose support is an iterated function system (IFS).
Constructor requires only an IFS, which is of type ```Array{Similarity}```.
All other essential properties can be deduced from this, including barycentre, diameter
and dimension, which are approximated numerically.

Has the outer constructor, which only requires IFS (a vector of [`Similarity`](@ref)) as an input.

    InvariantMeasure(sims::Vector{Similarity}; measure::Real=1.0) = 
        InvariantMeasure(sims, get_diameter(sims); measure=measure)
    
# Fields
- IFS: The iterated function system, a vector of [`Similarity`](@ref), describing the fractal
- spatial_dimension: The integer dimension of the smallest open set containing the fractal
- Hausdorff_dimension: The Hausdorff dimenson of the fractal
- homogeneous: A flag, true if all the contractions are of the same size
- Hausdorff_weights: A flag, true if this is a Hausdorff measure
- barycentre: The barycentre of the measure
- diameter: The diemater of the fractal support
- measure: The measure of the whole fractal, usually set to one
- weights: The probability weights describing the invariant measure
- disjoint: Flag for if the fractal support is disjoint
- connectedness: Matrix describing which subcomponents are connected
- symmetry_group: Vector of ```AutomorphicMap```, describing symmetries of the measure
"""
struct InvariantMeasure{V,M} <: SelfSimilarFractal{V,M}
    IFS::Vector{Similarity{V,M}}
    spatial_dimension::Int64
    Hausdorff_dimension::Float64
    homogeneous::Bool
    Hausdorff_weights::Bool
    barycentre::V
    diameter::Float64
    measure::Float64
    weights::Vector{Float64}
    disjoint::Bool
    connectedness::Matrix{Bool}
    symmetry_group::Vector{AutomorphicMap{V,M}}
end

# InvariantMeasure(sims::Array{Similarity}; measure::Real=1.0) = InvariantMeasure(sims, get_diameter(sims); measure=measure)

function InvariantMeasure(sims::Vector{Similarity{V,M_}}; diameter::Real=0.0, measure::Real=1.0,
        weights::Vector{<:Real}=[0.0], connectedness::Matrix{Bool}=Matrix(I(length(sims))),
        symmetry_group::Vector{AutomorphicMap{V,M_}}=TrivialGroup(length(sims[1].δ))) where {V<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}# outer constructor for parent_measure type
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
            weights[m] = sims[m].r^Hdim
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

    # disjointness
    connectedness==Matrix(I(length(sims))) ? disjoint=true : disjoint=false

    return InvariantMeasure(sims, top_dim, Hdim, homogeneous, Hausdorff_weights, Sbary,
            Float64(diameter), Float64(measure), weights, disjoint, connectedness, symmetry_group)

end

"""
    struct SubInvariantMeasure{V,M} <: SelfSimilarFractal{V,M}
        parent_measure::InvariantMeasure
        IFS::Vector{Similarity{V,M}} # could be removed ?
        index::Vector{Int64}
        barycentre::V
        diameter::Float64
        measure::Float64
    end

Represents a fractal measure which has been derived from an `InvariantMeasure`.

# Fields
- parent_measure: This measure is supported on a subet of the support of parent_measure
- IFS: The vector of similarities describing the fractal support of this measure
- index: The vector index corresponding to the contractions applied to parent_measure
- barycentre: The barycentre of this measure
- diameter: The diameter of the fractal support
- measure: The measure of the support
"""
struct SubInvariantMeasure{V,M} <: SelfSimilarFractal{V,M}
    parent_measure::InvariantMeasure
    IFS::Vector{Similarity{V,M}} # could be removed ?
    index::Vector{Int64}
    barycentre::V
    diameter::Float64
    measure::Float64
end
# outer constructor
"""
Representation of a subcomponent of a fractal Γ, using standard vector index notatation.
If Γ is a subattractor, then the vector indices are concatenated to produce a new subatractor,
which stores the original parent_measure.
"""
function SubInvariantMeasure(Γ::InvariantMeasure{V,M}, index::Vector{<:Integer}) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    #quick condition for trivial case:
    if index == [0]
        return SubInvariantMeasure{V,M}(Γ, Γ.IFS, [0], Γ.barycentre, Γ.diameter, Γ.measure)
    else #non-trivial case

        # get new measure and diameter. First initialise:
        new_diam = Γ.diameter
        new_measure = Γ.measure
        # now adjust the diameter and measure for the smaller scale:
        for m = index
            new_diam *= Γ.IFS[m].r
            new_measure *= Γ.weights[m]
        end

        new_IFS = Γ.IFS
        for m = index[end:-1:1]
            new_IFS = sim_map(Γ.IFS[m], new_IFS)
        end
        new_bary = get_barycentre(new_IFS, Γ.weights)

        # I think the below code is wrong... and I've been getting away with it somehow...
        # #start at old barycentre and map
        # new_bary = Γ.barycentre
        # for m=index[end:-1:1]
        #     new_bary = sim_map(Γ.IFS[m], new_bary)
        # end

        return SubInvariantMeasure{V,M}(Γ,new_IFS,   index, new_bary, new_diam, new_measure)
        
    end
end

function SubInvariantMeasure(Γ::SubInvariantMeasure{V,M}, index::Vector{<:Integer}) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
        #quick condition for trivial case:
        if index == [0]
            return Γ
        else #non-trivial case
    
            # get new measure and diameter. First initialise:
            new_diam = Γ.diameter
            new_measure = Γ.measure
            # now adjust the diameter and measure for the smaller scale:
            for m = index
                new_diam *= Γ.parent_measure.IFS[m].r
                new_measure *= Γ.parent_measure.weights[m]
            end
    
            #start as old barycentre and map
            # COULD BE IMPROVED, RE-APPLYING MAPS WHICH HAVE BEEN USED BEFORE
            new_IFS = Γ.IFS
            for m = index[end:-1:1]
                new_IFS = sim_map(Γ.IFS[m], new_IFS)
            end
            new_bary = get_barycentre(new_IFS, Γ.parent_measure.weights)


            # get new vector index and barycentre:
            if Γ.index != [0] # excluding trivial case
                index = vcat(Γ.index,index)
            end
    
    
            return SubInvariantMeasure{V,M}(Γ.parent_measure,new_IFS, index, new_bary, new_diam, new_measure)
        end
    end

SubInvariantMeasure(Γ::SelfSimilarFractal, index::Integer) = SubInvariantMeasure(Γ, [index])

# overload the indexing function, so we can get neat vector subscripts
getindex(Γ::SelfSimilarFractal, inds...) = SubInvariantMeasure(Γ,[i for i in inds])
getindex(Γ::SelfSimilarFractal, inds::Vector{<:Integer}) = SubInvariantMeasure(Γ,inds)

getweights(Γ::SelfSimilarFractal) = isa(Γ,InvariantMeasure) ? Γ.weights : Γ.parent_measure.weights
get_symmetry_group(Γ::SelfSimilarFractal) = isa(Γ,InvariantMeasure) ? Γ.symmetry_group : Γ.parent_measure.symmetry_group
get_spatial_dimension(Γ::SelfSimilarFractal) = isa(Γ,InvariantMeasure) ? Γ.spatial_dimension : Γ.parent_measure.spatial_dimension
get_connectedness(Γ::SelfSimilarFractal) = isa(Γ,InvariantMeasure) ? Γ.connectedness : Γ.parent_measure.connectedness
#get_barycentre(new_IFS, Γ.parent_measure.weights)
function changeweights(Γ::InvariantMeasure, μ::Vector{Float64}; G::Vector{AutomorphicMap{V,M}}=TrivialGroup(Γ.spatial_dimension)) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    return InvariantMeasure(Γ.IFS, Γ.spatial_dimension, Γ.Hausdorff_dimension, Γ.homogeneous, are_weights_Hausdorff(μ,Γ.IFS,Γ.Hausdorff_dimension), get_barycentre(Γ.IFS, μ), Γ.diameter, Γ.measure, μ, Γ.disjoint, Γ.connectedness, G)
end
function changeweights(Γ::SubInvariantMeasure, μ::Vector{Float64}; G::Vector{AutomorphicMap{V,M}}=TrivialGroup(length(sims[1].δ))) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    return changeweights(Γ.parent_measure,μ,G=G)[Γ.index]
end
#SubInvariantMeasure(changeweights(Γ.parent_measure,μ), Γ.IFS, Γ.index, Γ.barycentre, Γ.diameter, Γ.measure)

struct UnionInvariantMeasure{V,M} <: FractalMeasure{V,M}
    invariant_measures::Vector{InvariantMeasure{V,M}}
    spatial_dimension::Int64
end

function UnionInvariantMeasure(Γs::Vector{InvariantMeasure{V,M}}
            ) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    # first get max spatial dimension
    max_spatial_dim = 0
    for Γ ∈ Γs
        max_spatial_dim = max(max_spatial_dim,Γ.spatial_dimension)
    end

    # now embed each fractal in the max spatial dimension,
    # so everything lives in the same space
    for (n,Γ) = enumerate(Γs)
        if Γ.spatial_dimension<max_spatial_dim
            Γs[n] = embed(Γs[n], zeros(max_spatial_dim-Γ.spatial_dimension))
        end
    end

    # return the union, with consistent types
    return UnionInvariantMeasure(Γs, max_spatial_dim)
end