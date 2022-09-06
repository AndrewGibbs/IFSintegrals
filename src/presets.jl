# quick function which computes the barycentre, if it's not the default one with equal weights
function default_bary(S::Vector{Similarity{V,M}}, d::Number, weights::Vector{Float64}, Hausdorff_measure_bary::Union{Real,AbstractVector}) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    if are_weights_Hausdorff(weights,S,d)
        bary = V(Hausdorff_measure_bary)
    else
        bary = get_barycentre(S,weights)
    end
    return bary
end

"""
    CantorSet(;contraction = 1/3, weights=[1/2, 1/2])
Returns a Cantor Set as an Attractor (type) of an iterated function system.
"""
function CantorSet(;contraction = 1/3, weights=[1/2, 1/2])
    S = [Similarity(contraction,0.0),Similarity(contraction,1-contraction)]
    d = log(1/2)/log(contraction)
    bary = default_bary(S,d,weights,1/2)
    return Attractor(S,1,d,true,are_weights_Hausdorff(weights,S,d),bary,1.0,1.0,weights, true,Matrix(I(2)))
end

"""
    CantorDust(;contraction = 1/3, weights=[1/4, 1/4, 1/4, 1/4])
Returns a Cantor dust as an Attractor (type) of an iterated function system.
"""
function CantorDust(;contraction = 1/3, weights=[1/4, 1/4, 1/4, 1/4])
    S = [Similarity(contraction,[0.0,0.0]),Similarity(contraction,[1-contraction,0.0]),Similarity(contraction,[0.0,1-contraction]),Similarity(contraction,[1-contraction,1-contraction])]
    d = log(1/4)/log(contraction)
    bary = default_bary(S,d,weights,[0.5,0.5])
    return Attractor(S, 2, d, true, are_weights_Hausdorff(weights,S,d), bary, sqrt(2), 1.0,weights,true,Matrix(I(4)))
end


"""
    CantorN(N::Integer; contraction = 1/3)
Returns the Cartesian product of N Cantor Sets, as an Attractor (type) of an iterated function system.
For example, when N=2, we obtain Cantor Dust.
"""
function CantorN(N::Integer; contraction = 1/3)
    M = 2^N
    S = Similarity{SVector{N,Float64},SMatrix{N,N,Float64,N^2}}[]
    for m=1:M
        m_binary = digits(m-1,base=2,pad=N)
        push!(S,Similarity(contraction,m_binary.*(1-contraction)))
    end
    d = log(1/M)/log(contraction)
    weights = ones(M)./M 
    bary =  SVector{N,Float64}(0.5.*ones(N))
    return Attractor(S, 2, d, true, true, bary, sqrt(N), 1.0, weights, true, Matrix(Bool,I(2^N)))
end

"""
    Sierpinski(;weights=[1/3, 1/3, 1/3])
Returns the Sierpinski triangle, as an Attractor (type) of an iterated function system.
"""
function Sierpinski(;weights=[1/3, 1/3, 1/3])
    courage = Similarity(1/2,[0,0])
    wisdom = Similarity(1/2,[1/2,0])
    power = Similarity(1/2,[1/4,sqrt(3)/4])
    S = [courage,wisdom,power]
    d = log(3)/log(2)
    bary = default_bary(S,d,weights,[1/2,sqrt(3)/4])
    return Attractor(S, 2, d, true, are_weights_Hausdorff(weights,S,d), bary, 1.0, 1.0, weights,false,Matrix(ones(Bool,3,3)))
end
# function Sierpinski(;weights=[1/3, 1/3, 1/3])
#     courage = Similarity(1/2,[0,1/6])
#     wisdom = Similarity(1/2,sqrt(2)*[1/6,-1/6])
#     power = Similarity(1/2,sqrt(2)*[-1/6,-1/6])
#     S = [courage,wisdom,power]
#     d = log(3)/log(2)
#     bary = default_bary(S,d,weights,[0,(1-2*sqrt(2))/9])
#     return Attractor(S, 2, d, true, are_weights_Hausdorff(weights,S,d), bary, 1.0, 1.0, weights,false,Matrix(ones(Bool,3,3)))
# end
#            Attractor(sims,top_dim,Hdim,uniform,get_barycentre(sims,Hdim),diameter,measure)


"""
    SquareFlake(;weights=ones(16)./16)
Returns the Square Snowflake, sometimes referred to as the "Minkowski Island", as an Attractor (type) of an iterated function system.
See: https://en.wikipedia.org/wiki/Minkowski_sausage
"""
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
    return Attractor(IFS, 2, 2.0, true, are_weights_Hausdorff(weights,IFS,2), bary, R, 1.0, weights)
end

"""
    KochFlake(;weights = [1/3, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9])
Returns the Koch Snowflake as an Attractor (type) of an iterated function system.
"""
function KochFlake(;weights = [1/3, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9])
    IFS = [Similarity(sqrt(1/3),[0, 0], pi/6),
            Similarity(1/3,[0,2/3]),
            Similarity(1/3,[-1/sqrt(3),1/3]),
            Similarity(1/3,[-1/sqrt(3),-1/3]),
            Similarity(1/3,[0,-2/3]),
            Similarity(1/3,[1/sqrt(3),-1/3]),
            Similarity(1/3,[1/sqrt(3),1/3])
            ]
    bary = default_bary(IFS,2.0,weights,[0.0,0.0])
    # now create connectedness matrix for Γ_singularities
    M = 7
    connectendess = zeros(Bool,M,M)
    for m=1:M
        for m_=1:M
            if m==1 || m_==1
                connectendess[m,m_] = true
            end
            if abs(m-m_)==1 || abs(m-m_)==6
                connectendess[m,m_] = true
            end
        end
    end
    return Attractor(IFS, 2, 2.0, false, are_weights_Hausdorff(weights,IFS,2), bary, 2*sqrt(3)/3, 1.0, weights, false, connectendess)
end

function Carpet()
    ρ = 1/3
    IFS = [Similarity(ρ,[0,0]),
    Similarity(ρ,[0,ρ]),
    Similarity(ρ,[0,2ρ]),
    Similarity(ρ,[ρ,0]),
    Similarity(ρ,[ρ,2ρ]),
    Similarity(ρ,[2ρ,0]),
    Similarity(ρ,[2ρ,ρ]),
    Similarity(ρ,[2ρ,2ρ])
    ]
    M = 8
    connectendess = zeros(Bool,M,M)
    for m=1:M
        for m_=1:M
            if m==1 || m_==1
                connectendess[m,m_] = true
            end
            if abs(m-m_)==1
                connectendess[m,m_] = true
            end
        end
    end
    d = log(8)/log(3)
    bary = bary = default_bary(S,d,weights,[0.5,0.5])
    return Attractor(IFS, 2, d, true, are_weights_Hausdorff(weights,S,d), bary, sqrt(2), 1.0, weights, false, connectendess)
end