# General operations on rotations and reflections
struct GroupGenerator
    num_rotations :: Int64
    reflections :: Vector{Float64} # ∈ [0,π]
    centre :: Vector{Float64}  # ∈ R²
end

struct AutomorphicMap
    A::Matrix{Float64}
    δ::Vector{Float64}
end

# Now define some preset groups:
IdentityMap(n::Int64) = AutomorphicMap(Matrix(I(n)),zeros(n))
TrivialGroup(n::Int64) = [IdentityMap(n)]
# zero(AutomorphicMap{V}) = AutomorphicMap(Matrix(I(2)),[0.0,0.0])

# T₁(T₂x) = A₁(T₂x)+δ₁ = A₁(A₂x+δ₂)+δ₁ = A₁A₂x + A₁δ₂+δ₁
∘(T₁::AutomorphicMap, T₂::AutomorphicMap) = AutomorphicMap(T₁.A*T₂.A, T₁.A*T₁.δ+T₂.δ)

apply_automorphic_map(T::AutomorphicMap,x::Vector{Float64}) = T.A*x+T.δ

# stuff which constructs these group elements, using an angular argument
rotation2(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]
reflection2(θ) = [cos(2θ) sin(2θ); sin(2θ) -cos(2θ)]

# build full group, given the GroupGenerator object
function get_group_operations(G::GroupGenerator)
    δθ = 2π/G.num_rotations
    num_reflections = length(G.reflections)
    
    T = [IdentityMap(length(G.centre)) for _=1:(G.num_rotations+num_reflections)]

    counter = 0
    for θ=0:δθ:(2π-δθ)
        counter+=1
        Arot = rotation2(θ)
        T[counter] = AutomorphicMap(Arot, G.centre-Arot*G.centre)
    end
    for ϑ in G.reflections
        counter+=1
        Areflection = reflection2(ϑ)
        T[counter] = AutomorphicMap(Areflection, G.centre-Areflection*G.centre)
    end
    return T
end

# overload the approximation operator
isapprox(T₁::AutomorphicMap, T₂::AutomorphicMap) = T₁.A≈T₂.A && T₁.δ≈T₂.δ 

function DihedralGroup(n::Integer; centre = [0.0,0.0], angle_offest=0.0)
    if mod(n,2) == 0 # for even n, want to double frequency of reflections, but restrict to [0,π]
        δθ = π/n
        reflections = (0:δθ:(π-δθ)) .+ angle_offest
    else
        δθ = 2π/n # for odd n, split [0,2π] into n segments
        reflections = (0:δθ:(2π-δθ)) .+ angle_offest
    end
    return get_group_operations(GroupGenerator(n,reflections,centre))
end