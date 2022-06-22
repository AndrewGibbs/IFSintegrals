struct Projection{V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}} # onto the coefficient space of piecewise constants on the fractal
    domain::SelfSimilarFractal{V,M}
    #Lₕ::Vector{Vector{Int64}} # subindices list
    mesh::Vector{SubAttractor{V,M}}
    coeffs::Vector{<:Complex{<:Float64}}
end

function project(mesh::Vector{SubAttractor{V,M}}, f::Function, h_quad::Float64) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    F = zeros(Complex{Float64},length(mesh))
    n_count = 0
    for Γₙ in mesh
        n_count += 1
        x,w = barycentre_rule(Γₙ, h_quad)
        F[n_count] = w'*f.(x)
    end
    return Projection(mesh[1].attractor,mesh,F)
end

function \(K::DiscreteBIO, f::Function)
    fₕ = project(K.mesh, f, K.h_quad)
    coeffs = K.Galerkin_matrix \ fₕ.coeffs
    return Projection(K.BIO.domain, K.mesh, coeffs)
end

function \(K::DiscreteBIO, fₕ::Projection)
    thresh = 1E-8
    if length(fₕ.coeffs)>1000 #use iterative method
        # if K.BIO.self_adjoint && K.BIO.coercive
        #     coeffs = cg(K.Galerkin_matrix,fₕ.coeffs,reltol=thresh)
        # else
        if K.BIO.self_adjoint
            coeffs = minres(K.Galerkin_matrix,fₕ.coeffs,reltol=thresh)
        else
            coeffs = bicgstabl(K.Galerkin_matrix,fₕ.coeffs,1,reltol=thresh)
        end
    else
        coeffs = K.Galerkin_matrix \ fₕ.coeffs
    end
    return Projection(K.BIO.domain, K.mesh, coeffs)
end

# now some functions related to projections, and embeddings

function embed(f::Projection,g::Projection)
    #    make sure f is embedded into g
    if length(f.coeffs)>length(g.coeffs)
        F=f
        G=g
        f=G
        g=F
    end
    
    f_dim = length(f.coeffs)
    g_dim = length(g.coeffs)
    new_coeffs = zeros(Complex{Float64},g_dim)
    m_count = 1
    for mesh_el in g.mesh
        m = mesh_el.index
        # find the vector in f
        n_count = 0
        for nesh_el in f.mesh
            n = nesh_el.index
            n_count +=1
            if (n == m[1:length(n)]) || n == [0]
                new_coeffs[m_count] = f.coeffs[n_count]
                m_count+=1
                break
            end
        end
    end
    return Projection(f.domain,g.mesh,new_coeffs)
end

function embed(f::Projection,mesh::Vector{SubAttractor{V,M}}) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    if length(mesh) == length(f.coeffs)
        return f
    elseif length(mesh) > length(f.coeffs)
        f_dim = length(f.coeffs)
        mesh_dim = length(mesh)
        new_coeffs = zeros(Complex{Float64},mesh_dim)
        m_count = 1
        for mesh_el in mesh
            m = mesh_el.index
            # find the vector in f
            n_count = 0
            for nesh_el in f.mesh
                n = nesh_el.index
                n_count +=1
                if (n == m[1:length(n)]) || n == [0]
                    new_coeffs[m_count] = f.coeffs[n_count]
                    m_count+=1
                    break
                end
            end
        end
    else
        error("Mesh is not refined enough to support function")
    end
    return Projection(f.domain,mesh,new_coeffs)
end

function -(f::Projection,g::Projection)
    ϕ = embed(f,g)
    if length(f.coeffs)>length(g.coeffs)
        return Projection(ϕ.domain,ϕ.mesh,f.coeffs-ϕ.coeffs)
    else
        return Projection(ϕ.domain,ϕ.mesh,ϕ.coeffs-g.coeffs)
    end
end

function get_H_minus_half_norm_function(Γ::SelfSimilarFractal; h_BEM::Real=1.0, h_quad::Real=h_BEM, h_quad_diag::Real = h_quad,  vary_quad::Bool = true)
    Sᵢ = SingleLayer(Γ,im)
    DBIO = DiscreteBIO(Sᵢ; h_BEM = h_BEM, h_quad = h_quad, h_quad_diag = h_quad_diag, vary_quad = vary_quad)
    norm(ϕ::Projection) = sqrt((2*abs(embed(ϕ,DBIO.mesh).coeffs'*DBIO.Galerkin_matrix*embed(ϕ,DBIO.mesh).coeffs)))
    return norm
end

function get_H_minus_half_norm_function_from_matrix(Gᵢ::Matrix{ComplexF64},mesh::Vector{SubAttractor{V,M}}) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    #embed
    norm(ϕ::Projection) = sqrt((2*abs(embed(ϕ,mesh).coeffs'*Gᵢ*embed(ϕ,mesh).coeffs)))
    return norm
end

function get_H_minus_half_norm_function_from_matrix(Gᵢ::Matrix{ComplexF64})
    #embed
    norm(ϕ::Projection) = sqrt((2*abs(ϕ.coeffs'*Gᵢ*ϕ.coeffs)))
    return norm
end