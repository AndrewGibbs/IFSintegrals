using IFSintegrals

N = 20
h = 0.05

function rand_weights(N::Integer)
    w = rand(N)
    return w./sum(w)
end

ss = rand(10)
ss[1] = 0.0

Γs = [CantorSet(contraction = 1/2, weights = rand_weights(2)),
       CantorSet(weights = rand_weights(2)),
        KochFlake(),
        Sierpinski(),
        CantorDust()]

# GQ = (X,Y) -> gauss_quad(X,Y,N)
# CG = (X,Y) -> chaos_quad(X,Y,N)
function GQ(X::SelfSimilarFractal,Y::SelfSimilarFractal,f::Function)
    x,y,w = gauss_quad(X,Y,N)
    return w'*f.(x,y)
end

function CG(X::SelfSimilarFractal,Y::SelfSimilarFractal,f::Function)
    x,y,w = chaos_quad(X,Y,N)
    return w'*f.(x,y)
end

n_count = 0
for Γ ∈ Γs
    for s ∈ ss
        global n_count+=1
        println(n_count)
        Iₛ= s_energy(Γ,s,h)
        Iₛ = s_energy(Γ,s,h,μ₂ = rand_weights(length(Γ.IFS)))
        Iₛ = s_energy(Γ,s,CG)
        Iₛ = s_energy(Γ,s,CG,μ₂ = rand_weights(length(Γ.IFS)))
        if Γ.spatial_dimension == 1 # can also use Gauss quad rule here
            Iₛ = s_energy(Γ,s,GQ)
        end
    end
end