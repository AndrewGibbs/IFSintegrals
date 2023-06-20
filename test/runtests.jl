using IFSintegrals, LinearAlgebra, Test, MAT

# include("nb_egs_test.jl")
include("far_field_test.jl")
include("operator_tests.jl")
include("barycentre_rule_tests.jl")
include("NDJ_tests.jl")
include("gauss_tests.jl")
include("quad_bodge_test.jl")

screens = [CantorSet(), CantorDust(), Sierpinski(), KochFlake(),
            UnionInvariantMeasure([Sierpinski()+[2,2], CantorDust()])]
# UnionInvariantMeasure([CantorSet(),CantorSet()+1.5])
surfaces = [CantorDust(), UnionInvariantMeasure([CantorDust()+[2,2], CantorDust()])]
volumes = [KochFlake(), Dragon(), UnionInvariantMeasure([KochFlake()+[2,2], Dragon()])]
ss = rand(10)
h_energy = 0.1
operator_test_set = [SingleLayerOperatorLaplace(Sierpinski()),
                    SingleLayerOperatorLaplace(CantorDust(),ambient_dimension=3),
                    SingleLayerOperatorLaplace(CantorSet(),ambient_dimension=2)]

@testset "IFSintegrals tests" begin

    @testset "Lebesgue Far field comparison" begin
        @test get_cantor_far_field_err(8) ≈ 0 atol=1E-3
        @test get_cantor_far_field_err(11) ≈ 0 atol=1E-3
    end

    @testset "Monomial tests for gauss quad" begin
        @testset for Γ ∈ Γ_gauss, p = 1:50 # monomial_powers
        @test low_order_gauss(Γ,p) ≈ high_order_bary(Γ,p) rtol=1E-10
        end
    end

    @testset "surface BEMs" begin
        @testset for Γ ∈ surfaces
            @test(surface_test(Γ))
        end
    end

    @testset "screen BEMs" begin
        @testset  for Γ ∈ screens
            @test(screen_test(Γ))
        end
    end

    @testset "VIEs" begin
        @testset for Ω ∈ volumes
            @test(VIE_test(Ω))
        end
    end

    @testset "barycentre_rule_tests" begin
        @test(bary_test1())
        @test(bary_test2())
        @test(bary_test3())
        @test(bary_test4())
        # @test(bary_test5())
        @test(bary_test6())
        @test(bary_test7())
        @test(bary_test8())
    end

    @testset "Energy integral tests CG" begin
        @testset for s ∈ ss, Γ ∈ Γ_s_energy
                @test_nowarn(s_energy(Γ,s,h_energy))
                @test_nowarn(s_energy(Γ,s,h_energy,μ₂ = rand_weights(length(Γ.IFS))))
                @test_nowarn(s_energy(Γ, s, CG))
                @test_nowarn(s_energy(Γ, s, CG, μ₂ = rand_weights(length(Γ.IFS))))
            end
        end

    @testset "Energy integral test GQ" begin 
        @testset for s ∈ ss
        @test_nowarn(s_energy(CantorSet(),s,GQ))
        end
    end

    @testset "Laplace diagonal Galerkin entry near miss estimation" begin
        @testset for S ∈ operator_test_set, h_mesh ∈ [2.0, 0.5, 0.02]
                    @test get_first_matrix_entry(S,h_mesh,h_mesh/25)≈bodge_first_matrix_entry(S,h_mesh,h_mesh/100) rtol=2E-1
                end
    end

end;

# @testset "Testing all operators work with all fractals" begin
    #     @testset "surface BEMs" begin
    # @testset "surface BEMs" for Γ ∈ surfaces
    #     # @test_nowarn(surface_test(Γ))
    #     @test(true)
    # end
    #     end
    #     @testset "screen BEMs" begin
    #         @testset for Γ ∈ screens
    #             @test_nowarn(screen_test(Γ));
    #         end
    #     end
    #     @testset "VIEs" begin
    #         @testset for Ω ∈ volumes
    #             @test_nowarn(VIE_test(Ω));
    #         end
    #     end
    # end