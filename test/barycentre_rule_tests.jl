function get_test_fractal()
    ρ = 0.41 # uniform contraction
        IFS = [
            Similarity(ρ,[0,0])
            Similarity(ρ,[1-ρ,0])
            Similarity(ρ,[(1-ρ)/2,sqrt(3)*(1-ρ)/2])
            Similarity(ρ,[(1-ρ)/2,(1-ρ)/(2*sqrt(3))])
        ]
    return InvariantMeasure(IFS)
end


function test1(threshold = 1E-14)
    # don't change these parameters - these are specific to the benchmarks below
    h = 0.05
    Γ = get_test_fractal()
    y, w = barycentre_rule(Γ,h)

    f(x) = sin.(norm.(x)) # define integrand
    Ih = (w'*f(y))
    I_bm = 0.5593353808456712
    if abs(Ih-I_bm)>threshold
        return false
    else
        return true
    end
end

function test2(threshold = 1E-14)
    # don't change these parameters - these are specific to the benchmarks below
    h = 0.05
    Γ = get_test_fractal()
    Γ₁ = CantorSet(contraction = 1/3)
    Γ₂ = CantorSet(contraction = 3/8)
    f(x,y) = (sin.(abs.(x)).*cos.(abs.(y)))

    x,y,w = barycentre_rule(Γ₁,Γ₂,h)
    Ih = (w'*f(x,y))
    I_bm = 0.3727428300015248
    if abs(Ih-I_bm)>threshold
        return false
    else
        return true
    end
end

function test3(threshold = 1E-14)
    # don't change these parameters - these are specific to the benchmarks below
    h = 0.05
    Γ = get_test_fractal()
    t = 1/2 # singularity strength
    n = 2 #fixed point index, must be an integer
    Ih = eval_green_single_integral_fixed_point(Γ, t, h, n)
    I_bm = 1.4599772356543295
    if abs(Ih-I_bm)>threshold
        return false
    else
        return true
    end
end

function test4(threshold = 1E-14)
    # don't change these parameters - these are specific to the benchmarks below
    h = 0.05
    Γ = get_test_fractal()
    t = 1.0
    Ih = eval_green_double_integral(Γ, t, h)
    I_bm = 4.738268418342183
    if abs(Ih-I_bm)>threshold
        return false
    else
        return true
    end
end

function test5(threshold = 1E-14)
    # don't change these parameters - these are specific to the benchmarks below
    h = 0.05
    Γ = get_test_fractal()
    k = 1.96
    Ih = singular_elliptic_double_integral(Γ, k, h)
    I_bm = 0.3206113688339223 + 0.13651121576163053im
    if abs(Ih-I_bm)>threshold
        return false
    else
        return true
    end
end

function test6(threshold = 1E-14)
    # don't change these parameters - these are specific to the benchmarks below
    h = 0.05
    Γ = CantorSet()

    t = 0.0 # singularity strength
    n = 1 #fixed point index, must be an integer
    Ih = eval_green_single_integral_fixed_point(Γ, t, h, n)
    I_bm = -1.290945787479564
    if abs(Ih-I_bm)>threshold
        return false
    else
        return true
    end
end

function test7(threshold = 1E-14)
    # don't change these parameters - these are specific to the benchmarks below
    h = 0.05

    Γ = CantorSet()
    γ = Γ[1]
    t = 0.0 # singularity strength
    n = 1 #fixed point index, must be an integer
    Ih = eval_green_single_integral_fixed_point(γ, t, h, n)
    I_bm = -1.194257492236754
    if abs(Ih-I_bm)>threshold
        return false
    else
        return true
    end
end

function test8(threshold = 1E-14)
    # don't change these parameters - these are specific to the benchmarks below
    h = 0.05
    Γ = CantorDust()
    γ = Γ[1]
    t = 0.0 # singularity strength
    n = 1 #fixed point index, must be an integer
    Ih = eval_green_single_integral_fixed_point(γ, t, h, n)
    I_bm = -0.37960691168218147
    if abs(Ih-I_bm)>threshold
        return false
    else
        return true
    end
end