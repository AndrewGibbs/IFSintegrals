struct Projection{V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}} # onto the coefficient space of piecewise constants on the fractal
    domain::SelfSimilarFractal{V,M}
    Lₕ::Vector{Vector{Int64}} # subindices list
    # mesh::Vector{SubAttractor{V,M}}
    coeffs::Vector{<:Complex{<:Float64}}
end

function project(Γ::SelfSimilarFractal,Lₕ::Vector{Vector{Int64}},f::Function,h_quad::Float64)
    F = zeros(Complex{Float64},length(Lₕ))
    n_count = 0
    for n in Lₕ
        n_count += 1
        Γₙ = SubAttractor(Γ, n) # has already been computed elsewhere
        x,w = barycentre_rule(Γₙ, h_quad)
        F[n_count] = w'*f.(x)
    end
    return Projection(Γ,Lₕ,F)
end

function \(K::DiscreteBIO, f::Function)
    fₕ = project(K.BIO.domain, K.Lₕ, f, K.h_quad)
    coeffs = K.Galerkin_matrix \ fₕ.coeffs
    return Projection(K.BIO.domain, K.Lₕ, coeffs)
end

function \(K::DiscreteBIO, fₕ::Projection)
    thresh = 1E-8
    if length(fₕ.coeffs)>1000 #use iterative method
        if K.BIO.self_adjoint && K.BIO.coercive
            coeffs = cg(K.Galerkin_matrix,fₕ.coeffs,reltol=thresh)
        elseif K.BIO.self_adjoint
            coeffs = minres(K.Galerkin_matrix,fₕ.coeffs,reltol=thresh)
        else
            coeffs = bicgstabl(K.Galerkin_matrix,fₕ.coeffs,1,reltol=thresh)
        end
    else
        coeffs = K.Galerkin_matrix \ fₕ.coeffs
    end
    return Projection(K.BIO.domain, K.Lₕ, coeffs)
end

# now some functions related to projections, and embeddings

function embed(f,g)
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
    for m in g.Lₕ
        # find the vector in f
        n_count = 0
        for n in f.Lₕ
            n_count +=1
            if n == m[1:length(n)]
                new_coeffs[m_count] = f.coeffs[n_count]
                m_count+=1
                break
            end
        end
    end
    return Projection(f.domain,g.Lₕ,new_coeffs)
end

function -(f::Projection,g::Projection)
    ϕ = embed(f,g)
    if length(f.coeffs)>length(g.coeffs)
        return Projection(ϕ.domain,ϕ.Lₕ,f.coeffs-ϕ.coeffs)
    else
        return Projection(ϕ.domain,ϕ.Lₕ,ϕ.coeffs-g.coeffs)
    end
end