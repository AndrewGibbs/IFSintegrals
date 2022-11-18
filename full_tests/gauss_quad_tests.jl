using IFSintegrals

h = 10^-10

prob_weights = rand(2)
prob_weights = prob_weights./sum(prob_weights) # random weights for Cantor set

Γs = [CantorSet(),CantorSet(weights = prob_weights)]

scat_count = 0
for Γ ∈ Γs
    global scat_count = scat_count + 1
    println(scat_count)
    x_b,w_b = barycentre_rule(Γ,h) # prep high order barycentre rule

    for N=1:50
        p = 2N-1 # exactness of polynomial
        I_bary = w_b'*(x_b.^p) # compute barycentre approx
        x_g,w_g = gauss_quad(Γ,N) # get gauss rule
        I_gauss = w_g'*(x_g.^p) # compute lower order gauss approx
        err = abs(I_bary-I_gauss)/abs(I_bary)
        println(err)
    end

end