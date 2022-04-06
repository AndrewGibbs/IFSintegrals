function single_layer_potential(Γ::SelfSimilarFractal,k::Real,density::Projection,x::Vector{Float64}; h::Float64=0.1/k)
    # kernel is the same as the single layer BIO, so use that:
    Φₓ(Y::Vector{Float64}) = SingleLayer(Γ, k).kernel(x,[Y...,0.0])
    # have embedded n-dimensional y in n+1 dimensional Y
    val = 0.0
    for m_count = 1:length(density.Lₕ)
        m = density.Lₕ[m_count]
        Γₘ = SubAttractor(Γ,m)
        y,w = barycentre_rule(Γₘ,h)
        # now convert y to n+1 dimensions
        val += (w'*Φₓ.(y))*density.coeffs[m_count]
    end
    return val
end
# simplify for Laplace kernel case:
single_layer_potential(Γ::SelfSimilarFractal,density::Projection,points::Array{<:Real,2}) = single_layer_potential(Γ,0.0,density,points)

#have not tested far-field pattern yet:
function far_field_pattern(Γ::SelfSimilarFractal,k::Real,density::Projection,points::Array{<:Real,2};h=0.1/k)
    num_points = size(points)[1]
    vals = zeros(Complex{Float64},num_points)
    # far-field kernel:
    K(θ::Float64, ψ::Float64, x_1::Array{<:Real,1}, x_2::Array{<:Real,1}) = exp.(-im*k*(cos(θ)*x_1.+sin(ψ)*x_2))
    for m_count =1:length(density.Lₕ)
        m = density.Lₕ[m_count]
        Γₘ = SubAttractor(Γ,m)
        y,w = barycentre_rule(Γₘ,h)
        # now convert y to n+1 dimensions
        for n = 1:num_points
            vals[n] += w'*K(points[n,1],points[n,2],y[:,1],y[:,2])*density.coeffs[m_count]
        end
    end
    return vals
end

function far_field_pattern(Γ::SelfSimilarFractal,k::Real,density::Projection,points::Array{<:Real,1};h=0.1/k)
    num_points = length(points)
    vals = zeros(Complex{Float64},num_points)
    # far-field kernel:
    K(θ::Float64, x::Array{Float64,1}) = exp.(-im*k*(cos(θ)*x))
    for m_count =1:length(density.Lₕ)
        m = density.Lₕ[m_count]
        Γₘ = SubAttractor(Γ,m)
        y,w = barycentre_rule(Γₘ,h)
        # now convert y to n+1 dimensions
        for n = 1:num_points
            vals[n] += w'*K(points[n],y)*density.coeffs[m_count]
        end
    end
    return vals
end