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



