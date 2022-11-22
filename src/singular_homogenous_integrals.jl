"""
    eval_green_double_integral(Γ::Union{InvariantMeasure,SubAttractor}, t::Float64, h::Float64; μ::Vector{Float64} = Γ.weights)

Approximates the integral ∫_Γ∫_Γ Φ_t(x,y) dμ₁(y)dμ₂(x), where Φ_t is the Green's function for
the n-dimensional Laplace problem, and the integrals are with respect to Hausdorff measure.

The optional argument μ₂ will default to μ₁, which is the probability weights assigned to the measure of Γ.
But this can be specified manually, if it is required that both integrals are over a different measure.
"""
function eval_green_double_integral(Γ::SelfSimilarFractal, t::Real, h::Real; μ₂::Vector{Float64} = getweights(Γ))

    μ₁ = getweights(Γ)
    M = length(Γ.IFS)
    for μ = [μ₁,μ₂]
        length(μ)!=M ? error("μ must be a vector containing the same length as the IFS") : nothing
        sum(μ) ≈ 1 ? nothing : error("μ must sum to one")
    end

    log_sum = 0.0
    scale = 1.0
    smooth_integrals = 0.0

    if μ₁ == μ₂
        Γ_μ₂ = Γ
    else
        Γ_μ₂ = changeweights(Γ,μ₂)
    end

    for m=1:M
        scale -= Γ.IFS[m].r^(-t)*μ₁[m]*μ₂[m]
        if t == 0.0
            log_sum += Γ.measure^2 * μ₁[m]*μ₂[m] * log(Γ.IFS[m].r)
        end
        Γm = Γ[m]

        # can exploit symmetry of double sum

        if μ₁ == μ₂
            for n=(m+1):M
                if m!=n
                    X1, X2, W = barycentre_rule(Γm,Γ[n],h)
                    smooth_integrals += 2*W'* Φₜ.(t,X1,X2)
                end
            end
        else
            for n=1:M
                if m!=n
                    X1, X2, W = barycentre_rule(Γm,Γ_μ₂[n],h)
                    smooth_integrals += W'* Φₜ.(t,X1,X2)
                end
            end
        end
    end
    return (smooth_integrals + log_sum)/scale
end

"""
    eval_green_single_integral_fixed_point(Γ::Union{InvariantMeasure,SubAttractor}, t::Float64, h::Float64, n::Int64)
    
Approximates the integral ∫_Γ Φ_t(x,y) dμ(x), where Φ_t is the Green's function for
the n-dimensional Laplace equation, and the integrals are with respect to Hausdorff measure.
"""
function eval_green_single_integral_fixed_point(Γ::SelfSimilarFractal, t::Real, h::Real, n::Int64)
    if isa(Γ,InvariantMeasure)
        μ = Γ.weights
    else
        μ = Γ.attractor.weights
    end
    ηₙ = fixed_point(Γ.IFS[n]) # get fixed point, convert to standard vector type
    Φₜ_(x) = Φₜ(t,x,ηₙ)
    M = length(Γ.IFS)
    if t == 0.0
        log_sum = Γ.measure*μ[n]*log(Γ.IFS[n].r)
    else
        log_sum = 0.0
    end
    smooth_integrals = 0.0
    scale = 1.0 - (Γ.IFS[n].r^-t*μ[n])
    for m=1:M
        if m!=n
            x, w = barycentre_rule(Γ[m],h)
            smooth_integrals += w'* Φₜ_.(x)
        end
    end

    return (smooth_integrals + log_sum)/scale

end