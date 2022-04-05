"""
    x,w = barycentre_rule(Γ::Union{Attractor,SubAttractor},h::Real) 

returns N weights w ∈ Rⁿ and nodes x ∈ Rᴺˣⁿ,
for approximation of integrals defined on an IFS Γ
"""
# function barycentre_rule(Γ::SelfSimilarFractal,h::Real) #where {V<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
#     if Γ.homogeneous && Γ.Hausdorff_weights # uniform attractor, can do things a bit quicker
#         h = min(h,Γ.diameter)
#         ℓ = Int64(ceil(log(h/Γ.diameter)/log(Γ.IFS[1].r)))
#         x,w = barycentre_uniform(Γ,ℓ)
#     else
#         N = length(Γ.barycentre)
#         Lₕ = subdivide_indices(Γ,h)
#         if N==1
#             x = zeros(Float64, length(Lₕ))
#         else
#             x = [zeros(Float64, N) for _=1:length(Lₕ)]
#         end
#         w = zeros(Float64, length(Lₕ))
#         for j =1:length(Lₕ)
#             γ = SubAttractor(Γ,Lₕ[j])
#             if N==1
#                 x[j] = γ.barycentre[1]
#             else
#                 x[j] = γ.barycentre
#             end
#             w[j] = γ.measure
#         end
#     end

#     return x,w
# end

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
        for i=0:(n1-1)
            for j=0:(n2-1)
                X1[j*n1 + i + 1,:] = x1[i + 1,:]
                X2[j*n1 + i + 1,:] = x2[j + 1,:]
                W[j*n1 + i + 1] = w1[i + 1]*w2[j + 1]
            end
        end
    return X1, X2, W
end

# function barycentre_uniform(Γ::SelfSimilarFractal,ℓ::Int64)
#     # initialise arrays
#     N = length(Γ.barycentre)
#     M = length(Γ.IFS)
#     z = zeros(Float64, ℓ+1, N)
#     x = zeros(Float64, M^ℓ, N)
#     w = Γ.measure*Γ.IFS[1].r^(Γ.Hausdorff_dimension*ℓ)*ones(Float64,M^ℓ)
#     count = 0

#     if isa(Γ,Attractor)
#         z[1,:] = Γ.barycentre
#     else
#         z[1,:] = Γ.attractor.barycentre
#     end

#     m_inds = ones(Int64,ℓ)

#     for m = 1:M^ℓ
#         for j=1:ℓ
#             if mod(count,M^(ℓ-j))==0
#                 z[j+1,:] = sim_map(Γ.IFS[m_inds[j]], z[j,:])
#                 m_inds[j] += 1
#                 m_inds = mod.(m_inds .- 1,M) .+ 1
#             end
#         end
#         count += 1
#         x[count, :] = z[ℓ+1,:]
#     end
#     if isa(Γ,SubAttractor)
#         if Γ.index != [0]
#             for m = Γ.index[end:-1:1]
#                 x = sim_map(Γ.IFS[m], x)
#             end 
#             # X = slicematrix(x)
#         end
#     end
#     if N==1
#         return vec(x), w
#     else
#         return slicematrix(x), w
#     end
# end

function subdivide_indices(Γ::SelfSimilarFractal, h::Real; int_type::DataType=UInt64)
    I = Vector{int_type}[]
    M = length(Γ.IFS)
    r = zeros(M)

    if Γ.diameter >= h
        subdiv = true
        for m=int_type.(1:M)
            push!(I,[m])
            r[m] = Γ.IFS[m].r
        end
    else
        subdiv = false
    end

    while subdiv
        subdiv = false
        split_vecs = Int64[]
        for j = 1:length(I)
           if Γ.diameter*prod(r[I[j]]) >= h
                subdiv = true
                push!(split_vecs,j)
            end
        end
        if subdiv
            new_vecs = Vector{int_type}[]
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
        I = [[int_type(0)]]
    end
    return I#convert(Array{Array{Int64,1},1},I)
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

            # can exploit symmetry of double sum if attractor is uniform
            if true
                uniform_const = 2
                n_start = m+1
            else
                uniform_const = 1
                n_start = 1
            end

            for n=n_start:M
                Γn = SubAttractor(Γ,[n])
                if m!=n
                    X1, X2, W = barycentre_rule(Γm,Γn,h)
                    smooth_integrals += uniform_const*W'* Φₜ(t,X1,X2)
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
                smooth_integrals += w'* Φₜ(t,Matrix(x),ηₙ)
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

function chaos_quad(Γ::SelfSimilarFractal{V,M},N::Int64,x₀::AbstractVector) where {V<:AbstractVector, M<:Union{Real,AbstractMatrix}}
    x = [Vector{Float64}(undef, length(x₀)) for _ = 1:N]
    x[1] = x₀
    for n=2:N
        x[n] = chaos_step(Γ,x[n-1])
    end
    return x, ones(N)./N
end

function chaos_quad(Γ::SelfSimilarFractal,N::Int64,x₀::Float64)
    if Γ.topological_dimension>1
        error("Initial guess must match dimension of ambient space")
    end
    X,w = chaos_quad(Γ,N,[x₀])
    x = zeros(N)
    for n=1:N
        x[n] = X[n][1]
    end
    return x,w
end

function chaos_step(Γ::SelfSimilarFractal,x::Union{AbstractVector,Real})
    # randomly choose a map
    τ = minimum((1:length(Γ.weights))[rand().<cumsum(Γ.weights)])
    return sim_map(Γ.IFS[τ],x)
end