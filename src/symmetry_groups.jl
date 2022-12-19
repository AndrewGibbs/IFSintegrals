# General operations on rotations and reflections
struct GroupGenerator2D
    num_rotations :: Int64
    reflections :: Vector{Float64} # ∈ [0,π]
    centre :: Vector{Float64}  # ∈ R²
end

struct AutomorphicMap{V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    A::M
    δ::V
end

# Now define some preset groups:
function IdentityMap(ndims::Int64)
    if ndims == 1
        return AutomorphicMap(1.0,0.0)
    else
        return AutomorphicMap(SMatrix{ndims,ndims,Float64}(Matrix(1.0*I(ndims))), SVector{ndims,Float64}(zeros(ndims)))
    end
end

TrivialGroup(n::Int64) = [IdentityMap(n)]

# T₁(T₂x) = A₁(T₂x)+δ₁ = A₁(A₂x+δ₂)+δ₁ = A₁A₂x + A₁δ₂+δ₁
∘(T₁::AutomorphicMap, T₂::AutomorphicMap) = AutomorphicMap(T₁.A*T₂.A, T₁.A*T₁.δ+T₂.δ)

apply_automorphic_map(T::AutomorphicMap,x::Vector{Float64}) = T.A*x .+ T.δ

# stuff which constructs these group elements, using an angular argument
rotation2(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]
reflection2(θ) = [cos(2θ) sin(2θ); sin(2θ) -cos(2θ)]

# build full group, given the GroupGenerator object
function get_group_operations2D(G::GroupGenerator2D)
    δθ = 2π/G.num_rotations
    num_reflections = length(G.reflections)
    static_centre = SVector{2,Float64}(G.centre)
    
    T = [IdentityMap(2) for _=1:(G.num_rotations+num_reflections)]

    counter = 0
    for θ=0:δθ:(2π-δθ)
        counter+=1
        Arot = SMatrix{2,2,Float64}(rotation2(θ))
        T[counter] = AutomorphicMap(Arot, static_centre-Arot*static_centre)
    end
    for ϑ in G.reflections
        counter+=1
        Areflection = SMatrix{2,2}(reflection2(ϑ))
        T[counter] = AutomorphicMap(Areflection, static_centre-Areflection*static_centre)
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
    return get_group_operations2D(GroupGenerator2D(n,reflections,centre))
end

D₂_in_1D(;centre::Float64=0.0) = [IdentityMap(1), AutomorphicMap(-1.0, 2*centre)]