import LazySets: convex_hull, VPolygon, Singleton, element


function get_diam_long(sims::Array{Similarity}; tol=1E-8)
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

    X_slice = slicematrix(X')
    sX_slice = slicematrix(sX')
    # now determine how many iterations we need to reach prescribed tolerance
    num_its = Int(ceil(log(tol*(1-L)/(2h_dist(X_slice,sX_slice)))/log(L)))

    #quite a bit of converting back and forth between nested arrays and matrices here

    # now apply the full IFS map iteratively, until we are within diam tolerance
    Xₙ = convex_hull(sX_slice)
    for n=2:num_its
        X_ = hcat(Xₙ...)
        SX_ = full_map(sims,X_)
        SX_slice = slicematrix(SX_')
        Xₙ = convex_hull(SX_slice)
    end

    # diameter of hull at nth iteration is 
    return get_diameter(hcat(Xₙ...))
end

function get_diameter(sims::Array{Similarity})
    X = get_fixed_points(sims)
    X_slice = slicematrix(X')
    Xhull = VPolygon(convex_hull(X_slice))


    sX = full_map(sims,X)
    sX_slice = slicematrix(sX')

    X_in_hull = true
    for y in sX_slice
        element(Singleton(y)) ∈ Xhull ? true : X_in_hull = false
    end
    X_in_hull ? d=get_diameter(X) : d=get_diam_long(sims)
    return d
end

function get_diameter(X::Array{<:Real,2})
    M,_ = size(X)
    diam = 0
    for m=1:M
        for n=(m+1):M
            R = dist(X[m,:],X[n,:])
            if R>diam
                diam = R
            end
        end
    end
    return diam
end

function get_fixed_points(sims::Array{Similarity})
    M = length(sims)
    FPs = zeros(length(sims[1].δ),M)
    for m = 1:M
        FPs[:,m] = fixed_point(sims[m])
    end
    return FPs
end

"""Computes the Hausdorff distance between two sets"""
function h_dist(X::Array{<:Array{<:Real}},Y::Array{<:Array{<:Real}})
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