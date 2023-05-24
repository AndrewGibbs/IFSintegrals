using IFSintegrals, LinearAlgebra, Test

macro no_error(ex)
    quote
        try
            $(esc(ex))
            true
        catch
            false
        end
    end
end

function rand_weights(N::Integer)
    w = rand(N)
    return w./sum(w)
end
rand_contraction() = rand()/2

println("Barycentre rule tests")
include("barycentre_rule_tests.jl")
@testset begin
    @test(test1())
    @test(test2())
    @test(test3())
    @test(test4())
    @test(test5())
    @test(test6())
    @test(test7())
    @test(test8())
end

println("Energy integral tests")
include("NDJ_tests.jl")

ss = rand(10)
ss[1] = 0.0
N = 20
h = 0.05
Ncantordims=3

Γs = [CantorSet(contraction = 1/2, weights = rand_weights(2)),
       CantorSet(weights = rand_weights(2)),
        KochFlake(),
        Sierpinski(),
        CantorDust()]

@testset begin
    for Γ ∈ Γs
        for s ∈ ss
            # println(s,h)

            @test@no_error(s_energy(Γ,s,h))
            
            @test@no_error(s_energy(Γ,s,h,μ₂ = rand_weights(length(Γ.IFS))))
            
            @test@no_error(s_energy(Γ,s,CG))

            @test@no_error(s_energy(Γ,s,CG,μ₂ = rand_weights(length(Γ.IFS))))
            
            if Γ.spatial_dimension == 1
                @test@no_error(s_energy(Γ,s,GQ))
            end
        end
    end
end

println("Testing preset fractals")

Γs = [CantorSet(),
CantorSet(contraction = rand_contraction(), weights = rand_weights(2)),
CantorDust(),
CantorDust(contraction = rand_contraction(), weights = rand_weights(4)),
CantorN(Ncantordims),
CantorN(Ncantordims, contraction = rand_contraction()),
Sierpinski(),
Sierpinski(weights=rand_weights(3)),
SquareFlake(),
SquareFlake(weights=rand_weights(16)),
KochFlake(),
KochFlake(weights = rand_weights(7))]

q=2
Npts = 10

f(x) = sin(norm(x))

@testset begin
    for (n,Γ) ∈ enumerate(Γs)
        h = Γ.IFS[1].r^q*Γ.diameter
        H = [h-eps(), h, h+eps()]
        for h_ ∈ H
            @test@no_error(barycentre_rule(Γ,h))
            @test@no_error(barycentre_rule(Γ,Γ,h))
            @test@no_error(long_bary(Γ,f,h))
        end

        @test@no_error(chaos_quad(Γ, Npts))
        if Γ.spatial_dimension == 1
            @test@no_error(gauss_quad(Γ, Npts))
        end
    end
end

println("Testing Gauss quadrature")

h = 10^-8

prob_weights = rand(2)
prob_weights = prob_weights./sum(prob_weights) # random weights for Cantor set

Γs = [CantorSet(), CantorSet(weights = prob_weights)]

scat_count = 0

function exactness_check(Γ,N,x_b,w_b)
    p = 2N-1 # exactness of polynomial
    I_bary = w_b'*(x_b.^p) # compute barycentre approx
    x_g,w_g = gauss_quad(Γ,N) # get gauss rule
    I_gauss = w_g'*(x_g.^p) # compute lower order gauss approx
    # println(abs(I_bary-I_gauss)/abs(I_bary))
    close_enough = abs(I_bary-I_gauss)/abs(I_bary)<10^4*eps()
    return close_enough
end

@testset begin
    for Γ ∈ Γs
        # global scat_count = scat_count + 1
        # println(scat_count)
        x_b,w_b = barycentre_rule(Γ,h) # prep high order barycentre rule

        for N=1:50
            @test(exactness_check(Γ,N,x_b,w_b))
        end

    end
end

println("Testing notebook example files") 
include("nb_egs_test.jl")

println("Testing values of 2D Cantor set scattering")
include("far_field_test.jl")