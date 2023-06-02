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

Γ = Sierpinski()
Sₖ = SingleLayerOperatorLaplace(Γ,ambient_dimension=3)
h_quad = 0.005
Sₖₕ = DiscreteSIO(Sₖ, h_mesh=2.0, h_quad=h_quad)
current_code_val = Sₖₕ.Galerkin_matrix[1,1]
bodge_val = bary_bodge_singular_integral(Sₖₕ.mesh[1],Sₖₕ.SIO.kernel,h_quad)
# heatmap(abs.(Sₖₕ.Galerkin_matrix),yflip=true)