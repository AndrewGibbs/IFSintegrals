function get_quadrature_from_partition(γ::partition_data{T}, M::Int64, S::Vector{Similarity{T,M_}}, weights::Vector{Float64}, h::Float64) where {T<:Real, M_<:Real}
    X = create_partition(γ, M, S, weights, h)
    N = length(X)
    w = zeros(N)
    zero_node = zero(T)
    x = [zero_node for _=1:N]
    for n=1:N
        x[n] = X[n].barycentre
        w[n] = X[n].weight
    end
    return x, w
end

function get_quadrature_from_partition(γ::partition_data{T}, M::Int64, S::Vector{Similarity{T,M_}}, weights::Vector{Float64}, h::Float64) where {T<:AbstractVector, M_<:Union{Real,AbstractMatrix}}
    X = create_partition(γ, M, S, weights, h)
    N = length(X)
    w = zeros(N)
    zero_node = zeros(Float64,length(γ.barycentre))
    x = [zero_node for _=1:N]
    for n=1:N
        x[n] = Vector{Float64}(X[n].barycentre)
        w[n] = X[n].weight
    end
    return x, w
end

"""
    x,w = barycentre_rule(Γ::Union{Attractor,SubAttractor},h::Real) 

returns N weights w ∈ Rⁿ and nodes x ∈ Rᴺˣⁿ,
for approximation of integrals defined on an IFS Γ
"""
barycentre_rule(Γ::Attractor{T,M_},h::Float64) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}} = get_quadrature_from_partition(partition_data_unindexed{T}(Γ.barycentre,Γ.measure,Γ.diameter), length(Γ.IFS), Γ.IFS, Γ.weights, h)
barycentre_rule(Γ::SubAttractor{T,M_},h::Float64) where {T<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}} = get_quadrature_from_partition(partition_data_unindexed{T}(Γ.barycentre,Γ.measure,Γ.diameter), length(Γ.attractor.IFS), Γ.attractor.IFS, Γ.attractor.weights, h)


"""
    x,y,w = barycentre_rule(Γ₁::Union{Attractor,SubAttractor},Γ₂::Union{Attractor,SubAttractor},h::Real)

returns N weights w ∈ Rⁿ and nodes x,y ∈ Rᴺˣⁿ,
for approximation of double integrals over Γ₁,Γ₂.
"""
function barycentre_rule(Γ1::SelfSimilarFractal,Γ2::SelfSimilarFractal,h::Real)
        top_dims = Γ1.topological_dimension
        x1, w1 = barycentre_rule(Γ1,h)
        n1 = length(w1)
        x2, w2 = barycentre_rule(Γ2,h)
        n2 = length(w2)
        X1 = [zero(x1[1]) for _=1:n1*n2]
        X2 = [zero(x2[1]) for _=1:n1*n2]
        W = zeros(Float64, n1*n2)
        for i=0:(n1-1)
            for j=0:(n2-1)
                X1[j*n1 + i + 1] = x1[i + 1]
                X2[j*n1 + i + 1] = x2[j + 1]
                W[j*n1 + i + 1] = w1[i + 1]*w2[j + 1]
            end
        end
    return X1, X2, W
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