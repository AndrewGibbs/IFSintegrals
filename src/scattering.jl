function single_layer_potential(k::Real,density::IFSintegrals.Projection,h::Float64=0.1/k)
    W = Complex{Float64}[]
    Y = Vector{Float64}[]

    for m_count = 1:length(density.mesh)
        y,w = barycentre_rule(density.mesh[m_count],h)
        Y = vcat(Y,y)
        W = vcat(W,w*density.coeffs[m_count])
    end
    function R(x_::Vector{Float64})
        R_ = zeros(Float64,length(Y::Vector{Vector{Float64}}))
        for n=1:length((Y::Vector{Vector{Float64}}))
            R_[n] = norm([x_[1]-(Y::Vector{Vector{Float64}})[n][1],x_[2]-(Y::Vector{Vector{Float64}})[n][2],x_[3]])
        end
        return R_
    end
    Sϕ(x::Vector{Float64}) = (W::Vector{Complex{Float64}})'*(exp.(im*k*R(x))./(4π*(R(x))))
    return Sϕ
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