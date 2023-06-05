prob_weights = [1/3, 2/3]#rand(2)
prob_weights = prob_weights./sum(prob_weights) # random weights for Cantor set

Γ_gauss = [CantorSet(), CantorSet(weights = prob_weights)]

# function exactness_check(Γ,N,x_b,w_b)
#     p = 2N-1 # exactness of polynomial
#     I_bary = w_b'*(x_b.^p) # compute barycentre approx
#     x_g,w_g = gauss_quad(Γ,N) # get gauss rule
#     I_gauss = w_g'*(x_g.^p) # compute lower order gauss approx
#     # println(abs(I_bary-I_gauss)/abs(I_bary))
#     close_enough = abs(I_bary-I_gauss)/abs(I_bary)<10^4*eps()
#     return close_enough
# end

function high_order_bary(Γ::SelfSimilarFractal,p::Integer)
    h = 10^-9
    x_b,w_b = barycentre_rule(Γ,h)
    return w_b'*(x_b.^p)
end

function low_order_gauss(Γ::SelfSimilarFractal,p::Integer)
    N = ceil(Int64,(p+1)/2) # get min number of quad points for exactness
    x_g,w_g = gauss_quad(Γ,N) # get gauss rule
    return w_g'*(x_g.^p) # compute lower order gauss approx
end