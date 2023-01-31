function volume_potential(density::Projection{V,M}, k::Real; h_quad::Float64=0.1/k) where {V<:AbstractVector, M<:Union{Real,AbstractMatrix}}
    # get points on surface of Γ
    Y, W = full_surface_quad(density,h_quad)
    vec_length = length(W)

    dims = density.domain.spatial_dimension
    !(dims ≈  density.domain.Hausdorff_dimension) ? error("Hausdorff dimension must equal spatial dimension") : nothing
    # if dims == 2
        ϕ2(R::Float64) = HelhmoltzGreen2D(k,R)
    # elseif dims == 3
        ϕ3(R::Float64) = HelhmoltzGreen3D(k,R)
    # else
    #     error("Volume potential doesn't make sense for dimension not equal to 2 or 3.")
    # end

    function Nϕ(x::Vector{Float64})
        R = zeros(Float64,vec_length)
        for n=1:vec_length
            for j = 1:dims
                R[n] += (x[j]-Y[n][j]).^2
            end
        end
        R = Vector{Float64}(sqrt.(R))

        val = zero(ComplexF64)
        if dims == 2
            val = transpose(W::Vector{Complex{Float64}})*(ϕ2.(R))
        elseif dims == 3
            val = transpose(W::Vector{Complex{Float64}})*(ϕ3.(R))
        else
            error("Volume potential doesn't make sense for dimension not equal to 2 or 3.")
        end

        return val
    end
    return Nϕ
end