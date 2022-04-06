import LazySets: convex_hull, VPolygon, Singleton, element

"""
Most of the ideas in this section were adapted from results in
http://www.cse.dmu.ac.uk/~hamzaoui/papers/ifs.pdf

"""

function get_diam_long(sims::Array{Similarity{V,M}}; tol=1E-8) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    # get max Lipschitz constant of contraction
    L = 0
    for s in sims
        if s.r>L
            L=s.r
        end
    end

    # as first guess, take fixed points of Similarity
    X = get_fixed_points(sims)
    sX = full_map(sims,X)

    # now determine how many iterations we need to reach prescribed tolerance
    num_its = Int(ceil(log(tol*(1-L)/(2h_dist(X,sX)))/log(L)))

    #quite a bit of converting back and forth between nested arrays and matrices here

    # now apply the full IFS map iteratively, until we are within diam tolerance
    Xₙ = convex_hull(sX)
    for _=2:num_its
        SXₙ = full_map(sims,Xₙ)
        Xₙ = convex_hull(SXₙ)
    end
    # diameter of hull at nth iteration is 
    return get_diameter(Xₙ)
end

function get_diameter(sims::Vector{Similarity{V,M}}) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    FPs = get_fixed_points(sims)

    if V<:Real
        D = get_diameter(FPs)
    else
        Xhull = VPolygon(convex_hull(FPs))

        s_of_fixed_points = full_map(sims,FPs)

        X_in_hull = true
        for y in s_of_fixed_points
            element(Singleton(y)) ∈ Xhull ? true : X_in_hull = false
        end
        X_in_hull ? D=get_diameter(s_of_fixed_points) : D=get_diam_long(sims)
    end
    return D
end

get_diameter(X::Vector{<:Real}) = abs(maximum(X)-minimum(X))

function get_diameter(X::Vector{Vector{Float64}})
    M = length(X)
    # N = length(X[1])

    diam = 0
    for m=1:M
        for n=(m+1):M
            R = norm(X[m]-X[n])
            if R>diam
                diam = R
            end
        end
    end
    return diam
end

function get_fixed_points(sims::Array{Similarity{V,M_}}) where {V<:Real, M_<:Union{Real,AbstractMatrix}}
    M = length(sims)
    # define vector of vectors via list comprehension 
    FPs = [zero(V) for _ = 1:M]
    for m = 1:M
        FPs[m] = fixed_point(sims[m])
    end
    return FPs
end

function get_fixed_points(sims::Array{Similarity{V,M_}}) where {V<:AbstractVector, M_<:Union{Real,AbstractMatrix}}
    M = length(sims)
    N = length(sims[1].δ)
    # define vector of vectors via list comprehension 
    FPs = [Vector{Float64}(undef, N) for _ = 1:M]
    for m = 1:M
        FPs[m] = fixed_point(sims[m])
    end
    return FPs
end

"""Computes the Hausdorff distance between two sets"""
function h_dist(X::Vector{<:Vector{<:Real}},Y::Vector{<:Vector{<:Real}})
    D = zeros(length(X),length(Y))
    x_count = 0
    for x in X
        x_count += 1
        y_count = 0
        for y in Y
            y_count += 1
            D[x_count,y_count] = dist(x,y)
        end
    end

    max_x_minindex = 0
    max_x_mindist = 0.0
    for nx = 1:length(X)
        xY_mindist,xY_mindex = findmin(D[nx,:])
        if xY_mindist > max_x_mindist
            max_x_minindex = xY_mindex
            max_x_mindist = xY_mindist
        end
    end

    max_y_minindex = 0
    max_y_mindist = 0.0
    for ny = 1:length(Y)
        yX_mindist, yX_mindex = findmin(D[:,ny])
        if yX_mindist > max_y_mindist
            max_y_minindex = yX_mindex
            max_y_mindist = yX_mindist
        end
    end

    return max(max_x_mindist, max_y_mindist)
end