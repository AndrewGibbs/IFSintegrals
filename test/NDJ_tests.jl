# N = 20
# h = 0.05

# function rand_weights(N::Integer)
#     w = rand(N)
#     return w./sum(w)
# end

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

# @testset begin
#     for Γ ∈ Γs
#         for s ∈ ss
#             # println(s,h)
            
#             @check_it_runs(s_energy(Γ,s,h))
            
#             @check_it_runs(s_energy(Γ,s,h,μ₂ = rand_weights(length(Γ.IFS))))
            
#             @check_it_runs(s_energy(Γ,s,CG))

#             @check_it_runs(s_energy(Γ,s,CG,μ₂ = rand_weights(length(Γ.IFS))))
            
#             @check_it_runs(s_energy(Γ,s,GQ))
#         end
#     end
# end