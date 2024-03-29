# -------------- functions for affine operators on fractals --------------- #

# following function embeds similarity in equivalent one in higher/ambeint dimension
function embed(s::Similarity{V,M}, v::Vector{<:Real}
            ) where {V<:Union{AbstractVector}, M<:Union{Real,AbstractMatrix}}
    
    n = length(s.δ)
    n_ = n + length(v)
    A_ = Matrix{Float64}(I[1:n_,1:n_])
    if s.A != I
        A_[1:n,1:n] .= s.A
    end
    # return Similarity(s.r,vcat(s.δ,v),A_,rA_)
    return Similarity(s.r,vcat(s.δ,v),A_)
end

# function embed(s::Similarity{V,M}, v::Vector{<:Real}
#     ) where {V<:Union{AbstractVector}, M<:Union{Real,AbstractMatrix}}
    
#     n = length(s.δ)
#     n_ = n + length(v)
#     A_ = Matrix{Float64}(I[1:n_,1:n_])
#     rA_ = s.r*I[1:n_,1:n_]
#     if s.A != I
#     #     A_ = I
#     #     rA_ = s.r*I
#     # else
#         A_[1:n,1:n] .= s.A
#         rA_[1:n,1:n] .= s.rA
#     end
#     return Similarity(s.r,vcat(s.δ,v),A_,rA_)
# end

function embed(s::Similarity{V,M}, v::Vector{<:Real}
    ) where {V<:Real, M<:Real}
    if s.A == I
        A_ = I
        # rA_ = s.r*I
    else
        n = length(s.δ)
        n_ = n + length(v)
        A_ = Matrix{Float64}(I[1:n_,1:n_])
        A_[1:n,1:n] .= s.A
        # rA_ = s.r*I[1:n_,1:n_]
        # rA_[1:n,1:n] .= s.rA
    end
    # using Similarity constructor will make everything StaticArrays
    return Similarity(s.r,vcat(s.δ,v),A_)
end

function embed(AM::AutomorphicMap{V,M}, v::Vector{<:Real}
                ) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    n = length(AM.δ)
    n_ = n + length(v)
    A_ = Matrix{Float64}(I[1:n_,1:n_])
    A_[1:n,1:n] .= AM.A
    v_ = SVector{n_,Float64}(vcat(AM.δ,v))
    return AutomorphicMap(SMatrix{n_,n_,Float64,n_^2}(A_),v_)
end

function embed(Γ::InvariantMeasure, v::Vector{<:Real}
        )# where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}

    # get size of new vectors
    new_spat_dim = Γ.spatial_dimension + length(v)

    return InvariantMeasure{SVector{new_spat_dim,Float64},SMatrix{new_spat_dim,new_spat_dim,Float64,new_spat_dim^2}}(
        [embed(sₘ,v) for sₘ ∈ Γ.IFS], # embedded IFS
        new_spat_dim,
        Γ.Hausdorff_dimension,
        Γ.homogeneous,
        Γ.Hausdorff_weights,
        SVector{Γ.spatial_dimension + length(v), Float64}(vcat(Γ.barycentre,v)),
        Γ.diameter,
        Γ.measure,
        Γ.weights,
        Γ.disjoint,
        Γ.connectedness,
        [embed(T,v) for T ∈ Γ.symmetry_group]
    )
end
embed(Γ::InvariantMeasure, v::Real) = embed(Γ, [v])

# ------------------- addition / translation of fractals ----------------------

# following function translates fractal attractor
function +(s::Similarity, v::AbstractVector{<:Real})# where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    # if length(s.δ)<length(v)
    #     s = embed(s,zeros(length(v)))
    # end
    return Similarity(s.r, # scale (same)
                    -s.rA*v+s.δ+v, # translation 
                    s.A, # rotation/reflection (same)
                    s.rA) # scale*rotation/reflection (same)
end

+(AM::AutomorphicMap, v::AbstractVector{<:Real}) = AutomorphicMap(AM.A,-AM.A*v+AM.δ+v)

function +(Γ::InvariantMeasure{V,M}, v::AbstractVector{<:Real}
            ) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
            # make sure vector is right length, if not, mbed it in higher dimenson
            if Γ.spatial_dimension < length(v)
                Γ = embed(Γ,zeros(length(v)-Γ.spatial_dimension))
            elseif Γ.spatial_dimension > length(v)
                error("Spatial dimension of translation is lower than that of fractal.")
            end

            return InvariantMeasure(
                [sₘ+v for sₘ ∈ Γ.IFS], # translated IFS
                length(v),
                Γ.Hausdorff_dimension,
                Γ.homogeneous,
                Γ.Hausdorff_weights,
                SVector{length(v),Float64}(Γ.barycentre+v), # higher spatial dimension
                Γ.diameter,
                Γ.measure,
                Γ.weights,
                Γ.disjoint,
                Γ.connectedness,
                [T+v for T ∈ Γ.symmetry_group] # translated symmetry_group
            )
end

# sort out commutation & subtraction
+(v::AbstractVector{<:Real},s::Union{Similarity,AutomorphicMap, InvariantMeasure}) = +(s, v)
-(s::Union{Similarity,AutomorphicMap, InvariantMeasure}, v::AbstractVector{<:Real}) = +(s, -v)

# ------------- streching fractals ------------ #

function *(s::Similarity, ρ::Real)# where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    return Similarity(s.r, # scale (same)
                    ρ*s.δ, # translation 
                    s.A, # rotation/reflection (same)
                    s.rA) # scale*rotation/reflection (same)
end

*(AM::AutomorphicMap, ρ::Real) = AutomorphicMap(AM.A,AM.δ*ρ)

function *(Γ::InvariantMeasure{V,M}, ρ::Real,
    ) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    # make sure vector is right length, if not, mbed it in higher dimenson

    if ρ<=0
        error("scaling factor must be positive")
    end
    
    return InvariantMeasure(
        [sₘ*ρ for sₘ ∈ Γ.IFS], # translated IFS
        Γ.spatial_dimension,
        Γ.Hausdorff_dimension,
        Γ.homogeneous,
        Γ.Hausdorff_weights,
        Γ.barycentre*ρ,
        Γ.diameter*ρ,
        Γ.measure, # CONTENTIOUS ISSUE - leaving this unchanged
        Γ.weights,
        Γ.disjoint,
        Γ.connectedness,
        [T*ρ for T ∈ Γ.symmetry_group] # translated symmetry_group
    )
end

# communtes:
*(ρ::Real, Γ::InvariantMeasure) = *(Γ, ρ)

# ------------- rotating / reflecting fractals ------------ #

function *(T::Matrix{<:Real}, s::Similarity{V,M}
        ) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    if det(T) ≈ 1
        return Similarity(s.r, # scale (same)
            V(T*s.δ), # translation 
            M(T*s.A*inv(T)), # rotation/reflection
            M(T*s.rA*inv(T)) # scale*rotation/reflection
        )
    else
        error("rotation/reflection matrix needs to have determinant one")
    end
end

function *(T::Matrix{<:Real}, AM::AutomorphicMap{V,M}
        ) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    if det(T) ≈ 1
        n = length(AM.δ)
        return AutomorphicMap(
            M(T*AM.A*inv(T)), # rotation/reflection (same)
            V(T*AM.δ), # translation
        )
    else
        error("rotation/reflection matrix needs to have determinant one")
    end
end

function *(R::Matrix{<:Real}, Γ::InvariantMeasure{V,M},
    ) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    # make sure vector is right length, if not, mbed it in higher dimenson
    if !(det(R) ≈ 1)
        ρ = det(R)^(1/(size(R)[1]))
        Γ = ρ*Γ
        R = R/ρ
    end
    
    new_IFS = [R*sₘ for sₘ ∈ Γ.IFS] # translated IFS
    return InvariantMeasure(
        new_IFS,
        Γ.spatial_dimension,
        Γ.Hausdorff_dimension,
        Γ.homogeneous,
        Γ.Hausdorff_weights,
        SVector{length(Γ.barycentre),Float64}(R*Γ.barycentre), #get_barycentre(new_IFS,Γ.weights), # barycentre will change
        Γ.diameter,
        Γ.measure, # CONTENTIOUS ISSUE - leaving this unchanged
        Γ.weights,
        Γ.disjoint,
        Γ.connectedness,
        [R*T for T ∈ Γ.symmetry_group] # translated symmetry_group
    )
end