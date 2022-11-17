using IFSintegrals, LinearAlgebra

q = 3
Npts = 20
Ncantordims = 3
f(x) = sin(norm(x))

function rand_weights(N::Integer)
    w = rand(N)
    return w./sum(w)
end

rand_contraction() = rand()/2


 Γs = [CantorSet(),
    CantorSet(contraction = rand_contraction(), weights = rand_weights(2)),
    CantorDust(),
    CantorDust(contraction = rand_contraction(), weights = rand_weights(4)),
    CantorN(Ncantordims),
    CantorN(Ncantordims, contraction = rand_contraction()),
    Sierpinski(),
    Sierpinski(weights=rand_weights(3)),
    # SquareFlake(),
    # SquareFlake(weights=rand_weights(16)),
    KochFlake(),
    KochFlake(weights = rand_weights(7))]

# test quadrature rules

for (n,Γ) ∈ enumerate(Γs)
    println(n)
    h = Γ.IFS[1].r^q*Γ.diameter
    H = [h-eps(), h, h+eps()]
    for h_ ∈ H
        x,w = barycentre_rule(Γ,h)
        x,y,w = barycentre_rule(Γ,Γ,h)
        long_bary(Γ,f,h)
    end

    x,w = chaos_quad(Γ, Npts)
    if Γ.spatial_dimension == 1
        # x,w = gauss_quad(Γ, Npts)
    end
end