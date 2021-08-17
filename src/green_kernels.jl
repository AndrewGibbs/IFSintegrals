import SpecialFunctions: besselh

#pairwise distance between two arrays of N-dimensional points
function dist(x::Array{<:Real,2}, y::Array{<:Real,2})
    d1,d2 = size(x)
    #println(typeof(x))
    R = zeros(Float64,d1)
    for n=1:d2
        R += (x[:,n]-y[:,n]).^2
    end
    return sqrt.(R)
end

# distance between two single N-dimensional points
function dist(x::Array{<:Real,1}, y::Array{<:Real,1})
    d1 = length(x)
    R = 0.0
    for n=1:d1
        R += (x[n]-y[n])^2
    end
    return sqrt(R)
end

# Distances between each N-dimensional point of a vector, and a single N-dimensional point
function dist(x::Array{<:Real,1}, y::Array{<:Real,2})
    d1,d2 = size(y)
    #println(typeof(x))
    R = zeros(Float64,d1)
    for n=1:d2
        R += (x[n].-y[:,n]).^2
    end
    return sqrt.(R)
end
dist(x::Array{<:Real,2}, y::Array{<:Real,1}) = dist(y, x)
dist(x::Real,y::Real) = abs(x-y)

function Φₜ(t::Real, x::Union{Array{<:Real,2},Array{<:Real,1}}, y::Union{Array{<:Real,2},Array{<:Real,1}})
    if t==0
        return log.(dist(x,y))
    else
        return (dist(x,y)).^(-t)
    end
end

function zero_kernel(X::Union{Array{<:Real,2},Array{<:Real,1}}, Y::Union{Array{<:Real,2},Array{<:Real,1}})
    N = size(X)[1]
    return zeros(Complex{Float64},N)
end

HelhmoltzGreen2D(k::Real,x::Union{Array{<:Real,2},Array{<:Real,1}}, y::Union{Array{<:Real,2},Array{<:Real,1}}) = im/4*besselh.(0,1,k*dist(x,y))
function HelhmoltzGreen2D_Lipschitz_part(k::Real, X::Union{Array{<:Real,2},Array{<:Real,1}}, Y::Union{Array{<:Real,2},Array{<:Real,1}})
    N = size(X)[1]
    I = zeros(Complex{Float64},N)
    for n=1:N
        x = X[n,:]
        y = Y[n,:]
        if x == y
            I[n] = im/4 -1/(2π)*(0.577215664901532 + log(k/2))
        else
            I[n] = HelhmoltzGreen2D(k,x,y) + 1/(2π)*log(dist(x,y))
        end
    end
    return I
end

HelhmoltzGreen3D(k::Real,x::Union{Array{<:Real,2},Array{<:Real,1}},y::Union{Array{<:Real,2},Array{<:Real,1}}) = exp.(im*k*dist(x,y))./(4π*dist(x,y))
function HelhmoltzGreen3D_Lipschitz_part(k::Real,X::Union{Array{<:Real,2},Array{<:Real,1}}, Y::Union{Array{<:Real,2},Array{<:Real,1}})
    N = size(X)[1]
    I = zeros(Complex{Float64},N)
    for n=1:N
        x = X[n,:]
        y = Y[n,:]
        if x == y
            I[n] = im*k/(4π)
        else
            I[n] = HelhmoltzGreen3D(k,x,y) - 1.0 ./(4π*dist(x,y))
        end
    end
    return I
end