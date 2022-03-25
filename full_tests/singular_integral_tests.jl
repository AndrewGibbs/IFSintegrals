using IFSintegrals

ρ = 0.41 # uniform contraction
    IFS = [
        Similarity(ρ,[0,0])
        Similarity(ρ,[1-ρ,0])
        Similarity(ρ,[(1-ρ)/2,sqrt(3)*(1-ρ)/2])
        Similarity(ρ,[(1-ρ)/2,(1-ρ)/(2*sqrt(3))])
    ]
Γ = Attractor(IFS)


function test3(Γ,h)
    h = 0.05
    t = 1/2 # singularity strength
    n = 2 #fixed point index, must be an integer
    Ih = eval_green_single_integral_fixed_point(Γ, t, n)
    I_bm = 1.4599772356543295
    if abs(Ih-I_bm)>1E-12
        return false
    else
        return true
    end
end

function test4(Γ,h)
    k = 1.96
    Ih = singular_elliptic_double_integral(Γ, k, h)
    I_bm = 0.3206113688339223 + 0.13651121576163053im
    if abs(Ih-I_bm)>1E-12
        return false
    else
        return true
    end
end

print(test1())