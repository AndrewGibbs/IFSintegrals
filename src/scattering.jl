function full_surface_quad(density::Projection,h::Float64)
    W = Complex{Float64}[]
    Y = Vector{Float64}[]

    for m_count = 1:length(density.mesh)
        y,w = barycentre_rule(density.mesh[m_count],h)
        Y = vcat(Y,y)
        W = vcat(W,w*density.coeffs[m_count])
    end
    return Y,W
end

function single_layer_potential(k::Real,density::IFSintegrals.Projection, h_quad::Float64=0.1/k)
    Y, W = full_surface_quad(density,h_quad)
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

# far-field kernel 2D:
FF2D_kernel(θ::Float64, x::Float64,k::Float64) = exp(-im*k*(cos(θ)*x))

# far-field kernel 3D:
FF3D_kernel(θ::Float64, ψ::Float64, x::Vector{Float64}, k::Float64) = exp.(-im*k*(cos(θ)*x[1]+sin(ψ)*x[2]))

function far_field_pattern(ϕₕ::Projection, k::Real=0.0; h_quad::Float64=0.1/k)
    # I = 0.0
    # for (m,Γₘ) ∈ enumerate(ϕₕ.mesh)
    #     y,w = barycentre_rule(Γₘ,h_quad)
    #     I += w'*FF2D_kernel.(θ,y,k)*ϕₕ.coeffs[m]
    # end
    Y, W = full_surface_quad(ϕₕ,h_quad)
    if ϕₕ.domain.topological_dimension == 1
        Fϕ(θ::Float64) = transpose(W::Vector{Complex{Float64}})*FF2D_kernel.(θ,Y,k)
    elseif ϕₕ.domain.topological_dimension == 2
        Fϕ(θ::Float64, ψ::Float64) = transpose(W::Vector{Complex{Float64}})*FF3D_kernel.(θ,ψ,Y,k)
    end
    return Fϕ
end

# #have not tested far-field pattern yet:
# function far_field_pattern(Γ::SelfSimilarFractal,k::Real,density::Projection,points::Array{<:Real,2};h=0.1/k)
#     num_points = size(points)[1]
#     vals = zeros(Complex{Float64},num_points)
#     # far-field kernel:
#     K(θ::Float64, ψ::Float64, x_1::Array{<:Real,1}, x_2::Array{<:Real,1}) = exp.(-im*k*(cos(θ)*x_1.+sin(ψ)*x_2))
#     for m_count =1:length(density.Lₕ)
#         m = density.Lₕ[m_count]
#         Γₘ = SubAttractor(Γ,m)
#         y,w = barycentre_rule(Γₘ,h)
#         # now convert y to n+1 dimensions
#         for n = 1:num_points
#             vals[n] += w'*K(points[n,1],points[n,2],y[:,1],y[:,2])*density.coeffs[m_count]
#         end
#     end
#     return vals
# end

# function far_field_pattern(Γ::SelfSimilarFractal,k::Real,density::Projection,points::Array{<:Real,1};h=0.1/k)
#     num_points = length(points)
#     vals = zeros(Complex{Float64},num_points)
#     # far-field kernel:
#     K(θ::Float64, x::Array{Float64,1}) = exp.(-im*k*(cos(θ)*x))
#     for m_count =1:length(density.Lₕ)
#         m = density.Lₕ[m_count]
#         Γₘ = SubAttractor(Γ,m)
#         y,w = barycentre_rule(Γₘ,h)
#         # now convert y to n+1 dimensions
#         for n = 1:num_points
#             vals[n] += w'*K(points[n],y)*density.coeffs[m_count]
#         end
#     end
#     return vals
# end