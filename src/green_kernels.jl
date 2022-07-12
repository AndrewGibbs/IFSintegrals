import SpecialFunctions: besselh

function Φₜ(t::Real, x::T, y::T) where T<:Union{Real,AbstractVector}
    if t==0
        return log(norm(x-y))
    else
        return (norm(x-y))^(-t)
    end
end

Φₜ(t::Real, x::Vector{Float64}, y::SVector{N,Float64}) where N = Φₜ(t, x, Vector{Float64}(y))

zero_kernel(_, _) =  zero(Complex{Float64})

HelhmoltzGreen2D(k::Number,x::T, y::T) where T<:Union{Real,AbstractVector} = im/4*besselh(0,1,k*norm(x-y))
HelhmoltzGreen2D(k::Number,r::Real) = im/4*besselh(0,1,k*r)

function HelhmoltzGreen2D_Lipschitz_part(k::Number, x::T, y::T) where T<:Union{Real,AbstractVector}
    if x == y
        return im/4 -1/(2π)*(0.577215664901532 + log(k/2))
    else
        return HelhmoltzGreen2D(k,x,y) + 1/(2π)*log(norm(x-y))
    end
end

HelhmoltzGreen3D(k::Number,x::T,y::T) where T<:Union{Real,AbstractVector} = exp(im*k*norm(x-y))/(4π*norm(x-y))
HelhmoltzGreen3D(k::Number,r::Real) = exp(im*k*r)/(4π*r)

function HelhmoltzGreen3D_Lipschitz_part(k::Number, x::T, y::T) where T<:Union{Real,AbstractVector}
    if x == y
        return im*k/(4π)
    else
        return expm1(im*k*norm(x-y)) /(4π*norm(x-y))#HelhmoltzGreen3D(k,x,y) - 1.0 /(4π*norm(x-y))
    end
end

"""
    eval_green_double_integral(Γ::Union{Attractor,SubAttractor}, t::Float64, h::Float64; μ::Vector{Float64} = Γ.weights)

Approximates the integral ∫_Γ∫_Γ Φ_t(x,y) dμ(y)dμ'(x), where Φ_t is the Green's function for
the n-dimensional Laplace problem, and the integrals are with respect to Hausdorff measure.

The optional argument μ will default to μ', which is the probability weights assigned to the measure of Γ.
But this can be specified manually, if it is required that both integrals are over a different measure.
"""
function eval_green_double_integral(Γ::SelfSimilarFractal, t::Real, h::Real; μ₂::Vector{Float64} = [0.0])

    if isa(Γ,Attractor)
        μ₁ = Γ.weights
    else
        μ₁ = Γ.attractor.weights
    end
    if μ₂ == [0.0]
        μ₂ = μ₁
    end
    M = length(Γ.IFS)
    for μ = [μ₁,μ₂]
        length(μ)!=M ? error("μ must be a vector containing the same length as the IFS") : Nothing
        sum(μ) ≈ 1 ? Nothing : error("μ must sum to one")
    end

    log_sum = 0.0
    scale = 1.0
    smooth_integrals = 0.0

    for m=1:M
        scale -= Γ.IFS[m].r^(-t)*μ₁[m]*μ₂[m] #IFS[m].r^(2d-t)
        if t == 0.0
            log_sum += Γ.measure^2 * μ₁[m]*μ₂[m] * log(Γ.IFS[m].r) # Γ.measure^2 **IFS[m].r^(2d)*log(IFS[m].r)
        end
        Γm = Γ[m]#SubAttractor(Γ,[m])

        # can exploit symmetry of double sum

        for n=(m+1):M
            # Γn = SubAttractor(Γ,[n])
            if m!=n
                X1, X2, W = barycentre_rule(Γm,Γ[n],h)
                smooth_integrals += 2*W'* Φₜ.(t,X1,X2)
            end
        end
    end
    return (smooth_integrals + log_sum)/scale
end

"""
    eval_green_single_integral_fixed_point(Γ::Union{Attractor,SubAttractor}, t::Float64, h::Float64, n::Int64)
    
Approximates the integral ∫_Γ Φ_t(x,y) dH^d(x), where Φ_t is the Green's function for
the n-dimensional Laplace equation, and the integrals are with respect to Hausdorff measure.
"""
function eval_green_single_integral_fixed_point(Γ::SelfSimilarFractal, t::Real, h::Real, n::Int64)
    # if isa(Γ,Attractor)
    #     IFS = Γ.IFS
    #     d = Γ.Hausdorff_dimension
    # else
    #     IFS = Γ.attractor.IFS
    #     d = Γ.attractor.Hausdorff_dimension
    # end
    if isa(Γ,Attractor)
        μ = Γ.weights
    else
        μ = Γ.attractor.weights
    end
    ηₙ = Vector(fixed_point(Γ.IFS[n])) # get fixed point, convert to standard vector type
    Φₜ_(x) = Φₜ(t,x,ηₙ)
    M = length(Γ.IFS)
    if t == 0.0
        log_sum = Γ.measure*μ[n]*log(IFS[n].r) #Γ.measure*IFS[n].r^d*log(IFS[n].r)
    else
        log_sum = 0.0
    end
    smooth_integrals = 0.0
    # d = Γ.Hausdorff_dimension
    scale = 1.0 - (Γ.IFS[n].r^-t*μ[n])#1.0 - IFS[n].r^(d-t)
    for m=1:M
        # Γm = SubAttractor(Γ,[m])
        if m!=n
            x, w = barycentre_rule(Γ[m],h)
            smooth_integrals += w'* Φₜ_.(x)
        end
    end

    return (smooth_integrals + log_sum)/scale

end