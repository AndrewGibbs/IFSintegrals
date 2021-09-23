"""
    x,w = barycentre_rule(Γ::Union{Attractor,SubAttractor},h::Real) 

returns N weights w ∈ Rⁿ and nodes x ∈ Rᴺˣⁿ,
for approximation of integrals defined on an IFS Γ
"""
function barycentre_rule(Γ::Union{Attractor,SubAttractor},h::Real)
    Lₕ = subdivide_indices(Γ,h)
    # if Lₕ == [0] # this top condition may now be obsolete
    #     x = zeros(Float64, 1, Γ.topological_dimension)
    #     w = zeros(Float64, 1)
    #     x[1,:] = Γ.barycentre
    #     w = [Γ.measure]
    # else
        x = zeros(Float64, length(Lₕ), Γ.topological_dimension)
        w = zeros(Float64, length(Lₕ))
        for j =1:length(Lₕ)
            γ = SubAttractor(Γ,Lₕ[j])
            x[j,:] = γ.barycentre
            w[j] = γ.measure
        end
    # end
    return x,w
end

"""
    x,y,w = barycentre_rule(Γ₁::Union{Attractor,SubAttractor},Γ₂::Union{Attractor,SubAttractor},h::Real)

returns N weights w ∈ Rⁿ and nodes x,y ∈ Rᴺˣⁿ,
for approximation of double integrals over Γ₁,Γ₂.
"""
function barycentre_rule(Γ1::Union{Attractor,SubAttractor},Γ2::Union{Attractor,SubAttractor},h::Real)
        top_dims = Γ1.topological_dimension
        x1, w1 = barycentre_rule(Γ1,h)
        n1 = Int64(length(x1)/top_dims)
        x2, w2 = barycentre_rule(Γ2,h)
        n2 = Int64(length(x2)/top_dims)
        X1 = Array{Float64}(undef, n1*n2, top_dims)
        X2 = Array{Float64}(undef, n1*n2, top_dims)
        W = Array{Float64}(undef, n1*n2)
        # for j=0:(n1-1)
        #     for i=0:(n2-1)
        #         X1[j*n2 + i + 1,:] = x1[i + 1,:]
        #         X2[j*n2 + i + 1,:] = x2[j + 1,:]
        #         W[j*n2 + i + 1] = w1[i + 1]*w2[j + 1]
        #     end
        # end
        for i=0:(n1-1)
            for j=0:(n2-1)
                X1[j*n1 + i + 1,:] = x1[i + 1,:]
                X2[j*n1 + i + 1,:] = x2[j + 1,:]
                W[j*n1 + i + 1] = w1[i + 1]*w2[j + 1]
            end
        end
    return X1, X2, W
end

function subdivide_indices(Γ::Union{Attractor,SubAttractor},h::Real)
    I = []
    M = length(Γ.IFS)
    r = zeros(M)

    if Γ.diameter >= h
        subdiv = true
        for m=1:M
            push!(I,[m])
            r[m] = Γ.IFS[m].r
        end
    else
        subdiv = false
    end

    while subdiv
        subdiv = false
        split_vecs = []
        for j = 1:length(I)
           if Γ.diameter*prod(r[I[j]]) >= h
                subdiv = true
                push!(split_vecs,j)
            end
        end
        if subdiv
            new_vecs = []
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
        I = [[0]]
    end
    return convert(Array{Array{Int64,1},1},I)
end

"""
    eval_green_double_integral(Γ::Union{Attractor,SubAttractor}, t::Float64, h::Float64)

Approximates the integral ∫_Γ∫_Γ Φ_t(x,y) dH^d(y)dH^d(x), where Φ_t is the Green's function for
the n-dimensional Laplace problem, and the integrals are with respect to Hausdorff measure.
"""
function eval_green_double_integral(Γ::Union{Attractor,SubAttractor}, t::Real, h::Real; data=false)
    d = Γ.Hausdorff_dimension
    Npts = 0
    if t<d
        M = length(Γ.IFS)
        log_sum = 0.0
        scale = 1.0
        smooth_integrals = 0.0

        for m=1:M
        scale -= Γ.IFS[m].r^(2Γ.Hausdorff_dimension-t)
            if t == 0.0
                log_sum += Γ.measure^2*Γ.IFS[m].r^(2d)*log(Γ.IFS[m].r)
            end
            Γm = SubAttractor(Γ,[m])
            for n=1:M
                Γn = SubAttractor(Γ,[n])
                if m!=n
                    X1, X2, W = barycentre_rule(Γm,Γn,h)
                    smooth_integrals += W'* Φₜ(t,X1,X2)
                    Npts += length(W)
                end
            end
        end
        if data
            return Float64((smooth_integrals + log_sum)/scale), Npts
        else
            return Float64((smooth_integrals + log_sum)/scale)
        end
    else
        if data
            return Inf, 0
        else
            return Inf
        end
    end
end

"""
    eval_green_single_integral_fixed_point(Γ::Union{Attractor,SubAttractor}, t::Float64, h::Float64, n::Int64)
    
Approximates the integral ∫_Γ Φ_t(x,y) dH^d(x), where Φ_t is the Green's function for
the n-dimensional Laplace equation, and the integrals are with respect to Hausdorff measure.
"""
function eval_green_single_integral_fixed_point(Γ::Union{Attractor,SubAttractor}, t::Real, h::Real, n::Int64; data=false)
    d = Γ.Hausdorff_dimension
    Npts = 0
    if t < d
        M = length(Γ.IFS)
        if t == 0.0
            log_sum = Γ.measure*Γ.IFS[n].r^d*log(Γ.IFS[n].r)
        else
            log_sum = 0.0
        end
        smooth_integrals = 0.0
        d = Γ.Hausdorff_dimension
        ηₙ = fixed_point(Γ.IFS[n]) # get fixed point
        scale = 1.0 - Γ.IFS[n].r^(d-t)
        for m=1:M
            Γm = SubAttractor(Γ,[m])
            if m!=n
                x, w = barycentre_rule(Γm,h)
                smooth_integrals += w'* Φₜ(t,x,ηₙ)
                Npts += length(x)
            end
        end

        if data
            return Float64((smooth_integrals + log_sum)/scale), Npts
        else
            return Float64((smooth_integrals + log_sum)/scale)
        end
    else
        if data
            return Inf, 0
        else
            return Inf
        end
    end

end