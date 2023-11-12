# zero kernel, useful sometimes
zero_kernel(_, _) =  zero(Complex{Float64})

# Kernel for 2D Helmholtz BIE
HelhmoltzGreen2D(k::Number,r::Real) = im/4*besselh(0,1,k*r)
HelhmoltzGreen2D(k::Number,x::Real, y::Real) = im/4*besselh(0,1,k*abs(x-y))
HelhmoltzGreen2D(k::Number,x::AbstractVector, y::AbstractVector) = im/4*besselh(0,1,k*norm(x-y))

# Kernel for Laplace BIE, covers all dimensions, generalised from t=0,1 to any t
function Φₜ(t::Real, r::T) where T<:Union{Real,AbstractVector}
    if t==0
        return log(norm(r))
    else
        return norm(r)^(-t)
    end
end
Φₜ(t::Real, x::Real, y::Real)  = Φₜ(t,abs(x-y))
Φₜ(t::Real, x::AbstractVector, y::AbstractVector) = Φₜ(t,norm(x-y))

function HelhmoltzGreen2D_Lipschitz_part(k::Number, r::Real)# where T<:Union{Real,AbstractVector}
    if isapprox(r,0.0,atol=1e-14) # r ≈ 0
        return im/4 -1/(2π)*(0.577215664901532 + log(k/2))
    else
        return HelhmoltzGreen2D(k,r) + 1/(2π)*log(r)
    end
end
HelhmoltzGreen2D_Lipschitz_part(k::Number, x::Real, y::Real) = HelhmoltzGreen2D_Lipschitz_part(k,abs(x-y))
HelhmoltzGreen2D_Lipschitz_part(k::Number, x::AbstractVector, y::AbstractVector) = HelhmoltzGreen2D_Lipschitz_part(k,norm(x-y))

# 3D Helmholtz kernels
HelhmoltzGreen3D(k::Number,r::Real) = exp(im*k*r)/(4π*r)
HelhmoltzGreen3D(k::Number,x::Real,y::Real) = exp(im*k*abs(x-y))/(4π*abs(x-y))
HelhmoltzGreen3D(k::Number,x::AbstractVector,y::AbstractVector) = exp(im*k*norm(x-y))/(4π*norm(x-y))

function HelhmoltzGreen3D_Lipschitz_part(k::Number, r::Real)
    if isapprox(r,0.0,atol=1e-14)
        return im*k/(4π)
    else
        return expm1(im*k*norm(r)) /(4π*norm(r))
    end
end
HelhmoltzGreen3D_Lipschitz_part(k::Number, x::Real, y::Real) = HelhmoltzGreen3D_Lipschitz_part(k,abs(x-y))
HelhmoltzGreen3D_Lipschitz_part(k::Number, x::AbstractVector, y::AbstractVector) = HelhmoltzGreen3D_Lipschitz_part(k,norm(x-y))
