# idea here is to do a mismatched barycentre approx of a double integral
using IFSintegrals

function bary_bodge_singular_integral(Γ,f,h_quad)
    x1,w1 = barycentre_rule(Γ,h_quad)
    x2,w2 = barycentre_rule(Γ,h_quad*Γ.IFS[1].r*0.99999)

    n1 = length(w1)
    n2 = length(w2)
    X1 = repeat(x1,inner=n2)
    X2 = repeat(x2,outer=n1)
    W = repeat(w1,inner=n2).*repeat(w2,outer=n1)

    return W'*f.(X1,X2)
end

function get_first_matrix_entry(S,h_mesh,h_quad)
    Γ = S.domain
    Sₕ = DiscreteSIO(S, h_mesh=h_mesh, h_quad=h_quad)
    return Sₕ.Galerkin_matrix[1,1]
end

function bodge_first_matrix_entry(S,h_mesh,h_quad)
    Γ = S.domain
    mesh = IFSintegrals.mesh_fractal(Γ::InvariantMeasure, h_mesh::Number)
    return bary_bodge_singular_integral(mesh[1],S.kernel,h_quad)*S.singularity_scale
end

S₀ = SingleLayerOperatorLaplace(Sierpinski(),ambient_dimension=3)
h_mesh = 2.0
h_quad = h_mesh/25

≈(get_first_matrix_entry(S₀,h_mesh,h_quad), bodge_first_matrix_entry(S₀,h_mesh,h_quad), rtol=1e-1)
# heatmap(abs.(Sₖₕ.Galerkin_matrix),yflip=true)