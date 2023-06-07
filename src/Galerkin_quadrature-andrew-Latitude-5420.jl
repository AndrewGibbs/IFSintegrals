# The next couple of functions are designed to use less quadrature points in the BEM elements
# where this won't affect the accuracy.

F_nomeasure(r::Real, k::Number, n::Int64) = (1+(abs(k)*r)^(n/2+1))/r^(n+1)

function get_quad_scales(k::Real, spat_dim::Int64,mesh::Vector{<:SubInvariantMeasure})
    # compute upper and lower bounds for the F in my notes, which is stated above.
    N = length(mesh)
    F_upper = ones(Float64,N,N)
    F_lower = ones(Float64,N,N)
    for m_count = 1:N
        Γₘ = mesh[m_count]
        for n_count = m_count:N
            Γₙ = mesh[n_count]
            if n_count!=m_count
                dist_upper = norm(Γₙ.barycentre-Γₘ.barycentre)
                dist_lower = max(dist_upper - Γₙ.diameter - Γₘ.diameter,0)
                measure_weight = Γₙ.measure*Γₘ.measure
                # noting that F_nomeasure is monotonic decreasing in r, can bound as follows:
                if dist_lower>0
                    F_upper[m_count,n_count] = measure_weight*F_nomeasure(dist_lower, k, spat_dim)
                else
                    F_upper[m_count,n_count]= Inf
                end
                 F_lower[m_count,n_count] = measure_weight*F_nomeasure(dist_upper, k, spat_dim)
            else
                F_upper[m_count,n_count] = Inf
            end
            
            # now by symmetry
            F_upper[n_count,m_count] = F_upper[m_count,n_count]
            F_lower[n_count,m_count] = F_lower[m_count,n_count]
        end
    end
    # now can get a lower bound estimate on the quantity from my notes:
    quad_scales = Float16.(floor.(max.(sqrt.(maximum(F_lower)./F_upper),1),sigdigits = 4))
    # have made a change which rounds down to nearest Float16, to save space
    return quad_scales
end

function get_quad_scales(k::Real, spat_dim::Int64,
                        mesh1::Vector{<:SubInvariantMeasure},
                        mesh2::Vector{<:SubInvariantMeasure})
    N = length(mesh1)
    M = length(mesh2)
    F_upper = ones(Float64,M,N)
    F_lower = ones(Float64,M,N)

    for m = 1:M
        Γₘ = mesh1[m]
        for n = 1:N
            Γₙ = mesh2[n]
            dist_upper = norm(Γₙ.barycentre-Γₘ.barycentre)
            dist_lower = max(dist_upper - Γₙ.diameter - Γₘ.diameter,0)
            measure_weight = Γₙ.measure*Γₘ.measure
            # noting that F_nomeasure is monotonic decreasing in r, can bound as follows:
            if dist_lower>0
                F_upper[m,n] = measure_weight*F_nomeasure(dist_lower, k, spat_dim)
            else
                F_upper[m,n]= Inf
            end
            F_lower[m,n] = measure_weight*F_nomeasure(dist_upper, k, spat_dim)
        end
    end
    # now can get a lower bound estimate on the quantity from my notes:
    quad_scales = Float16.(floor.(max.(sqrt.(maximum(F_lower)./F_upper),1),sigdigits = 4))
    # have made a change which rounds down to nearest Float16, to save space
    return quad_scales
end