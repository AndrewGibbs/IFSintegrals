using IFSintegrals

function rand_weights(N::Integer)
    w = rand(N)
    return w./sum(w)
end

rand_contraction() = rand()/2

CantorSet()
CantorSet(contraction = rand_contraction(), weights = rand_weights(2))

CantorDust()
w = rand(4)
w = w./sum(w)
CantorSet(contraction = rand_contraction(), weights = rand_weights(4))

CantorN(4)
CantorN(4, contraction = rand_contraction())

Sierpinski()
Sierpinski(weights=rand_weights(3))

SquareFlake()
SquareFlake(weights=rand_weights(16))

KochFlake()
KochFlake(weights = rand_weights(7))