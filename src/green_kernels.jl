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