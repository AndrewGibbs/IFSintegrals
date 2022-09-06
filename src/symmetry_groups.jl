# General operations on rotations and reflections
struct GroupGenerator
    num_rotations :: Int64
    reflections :: Vector{Float64} # ∈ [0,π]
    centre :: Vector{Float64} # ∈ Rⁿ
end

struct AutomorphicMap
    A:Matrix{Float64}
    δ::Vector{Float64}
end

zero(AutomorphicMap) = AutomorphicMap(Matrix(I(2)),[0.0,0.0])

# T₁(T₂x) = A₁(T₂x)+δ₁ = A₁(A₂x+δ₂)+δ₁ = A₁A₂x + A₁δ₂+δ₁
∘(T₁::AutomorphicMap, T₂::AutomorphicMap) = AutomorphicMap(T₁.A*T₂*A, T₁.A*T₁.δ+T₂.δ)

apply_automorphic_map(T::AutomorphicMap,x::Vector{Float64}) = T.A*x+T.δ

# stuff which constructs these group elements, using an angular argument
rotation2(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]
reflection2(θ) = [cos(2θ) -sin(2θ); sin(2θ) - cos(2θ)]

# build full group, given the GroupGenerator object
function get_group_operations(G::GroupGenerator)
    δθ = 2π/G.num_rotations
    num_reflections = length(G.reflections)
    
    Trotation = zeros(AutomorphicMap,G.num_rotations)
    Treflection = zeros(AutomorphicMap,num_reflections)
    T = zeros(AutomorphicMap,num_rotations*num_reflections)

    for θ=0:δθ:(2π-δθ)
        Arot = rotation2(θ)
        Trotation = AutomorphicMap(Arot, G.centre-Arot*G.centre)
    end

    for ϑ in G.reflections
        Areflection = reflection2(ϑ)
        Treflection = AutomorphicMap(Arot, G.centre-Areflection*G.centre)
    end

    counter = 0
    for n_rot=1:G.num_rotations
        for n_refl=1:num_reflections
            counter+=1
            T[counter] = Trotation[counter]∘Treflection[counter]
        end
    end
    return T
end