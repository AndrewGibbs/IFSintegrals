# The below variation of vcat is needed for our convention that Γ₀:=Γ
vcat_(x,y) = vcat(x[x .!= 0],y[y .!= 0])
hcat_(x,y) = hcat(x[x .!= 0],y[y .!= 0])

function check_if_similar(Γ::SelfSimilarFractal{V,M},
    m::Vector{Int64},n::Vector{Int64},m_::Vector{Int64},n_::Vector{Int64},
    G::Vector{AutomorphicMap{V,M}}, G_::Vector{AutomorphicMap{V,M}}) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    # get shorthand for IFS
    S = Γ.IFS
    test_one = false
    test_two = false
    test_three = false
    is_similar = false # requires above three to be true
    ρ = 0.0 # initialise

    # define identity similarity, which is a workaround for index [0]
    s₀ = Similarity(1.0,zero(V),one(M),one(M))

    m != [0] ? sₘ = sim_comp(S,m) : sₘ = s₀
    n != [0] ? sₙ = sim_comp(S,n) : sₙ = s₀
    m_ !=[0] ? sₘ_ = sim_comp(S,m_) : sₘ_ = s₀
    n_ !=[0] ? sₙ_ = sim_comp(S,n_) : sₙ_ = s₀

    # first test (6)
    if isapprox(sₘ.r/sₙ.r, sₘ_.r/sₙ_.r, atol=100*eps())
        ρ = sₘ.r/sₙ.r
        test_one = true
    end

    if test_one
        for T ∈ G, T_ ∈ G_
            # second test (7)
            if isapprox(sₘ.A*T.A*inv(sₙ.A) , sₘ_.A*T_.A*inv(sₙ_.A), atol=100*eps())
                test_two = true
            else
                test_two = false
            end
            
            # third test (8)
            if isapprox(sₘ.δ .- sₘ_.δ .- sₘ.r*sₘ.A*(T.A/sₙ.r*inv(sₙ.A)*sₙ.δ .- T.δ), .- sₘ_.r*sₘ_.A*(T_.A/sₙ_.r*inv(sₙ_.A)*sₙ_.δ .- T_.δ), atol=100*eps())
                test_three = true
            else
                test_three = false
            end

            # if tests are satisfies for some particular group element, can break loop early:
            if test_two && test_three
                break
            end
        end
        is_similar = test_one && test_two && test_three
    end
    return is_similar, ρ
end

function check_for_similar_integrals(Γ::SelfSimilarFractal{V,M},
    X::Vector{Tuple{Vector{Int64}, Vector{Int64}}}, 
    mcat::Vector{Int64}, mcat_::Vector{Int64},
    G₁::Vector{AutomorphicMap{V,M}}, G₂::Vector{AutomorphicMap{V,M}},
    fubini_flag::Bool) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    is_X_similar = false
    similar_index = nothing
    proportionality_const = 0.0
    # should compactly write the following as a function, it's almost repeated
    for ∫∫_index = 1:length(X)
        ∫∫_indices_higher_level = X[∫∫_index]
        this_is_X_similar = true
        ρ = 1.0
        # else # check for less obvious similarities
        this_is_X_similar, ρ = check_if_similar(Γ, ∫∫_indices_higher_level[1], mcat, ∫∫_indices_higher_level[2], mcat_, G₁, G₂)
        # end
        if !this_is_X_similar && fubini_flag
            this_is_X_similar, ρ = check_if_similar(Γ, ∫∫_indices_higher_level[2], mcat, ∫∫_indices_higher_level[1], mcat_, G₁, G₂)
        end
        # if we've found a similarity, terminate the process early
        if this_is_X_similar
            is_X_similar = true
            proportionality_const = ρ
            similar_index = ∫∫_index
            break
        end
    end
    return is_X_similar, proportionality_const, similar_index
end

function check_for_ℓ_singular_integrals(Γ::SelfSimilarFractal{V,M_}, mcat::Vector{Int64}, mcat_::Vector{Int64}) where {V<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}
    is_singular = false

    # get important bits
    M = length(Γ.IFS)
    Γ_singularities = Γ.connectedness
    ℓ_depth = Int64(round(log(size(Γ_singularities)[1])/log(M)))

    if length(mcat) == length(mcat_) <= ℓ_depth
        ℓ = length(mcat)
        mentry = 0
        m_entry = 0
        for ℓ_=1:(ℓ-1)
            mentry += M^(ℓ-ℓ_)*(mcat[ℓ_]-1)
            m_entry += M^(ℓ-ℓ_)*(mcat_[ℓ_]-1)
        end
        mentry += mcat[end]
        m_entry += mcat_[end]
        Γ_singularities[mentry,m_entry] ? is_singular = true : nothing
    end

    return is_singular
end

function construct_singularity_matrix(Γ::SelfSimilarFractal{V,M_}, s::Number; μ₂::Vector{Float64} = getweights(Γ), G₂::Vector{AutomorphicMap{V,M_}}=get_symmetry_group(Γ)) where {V<:Union{Real,AbstractVector}, M_<:Union{Real,AbstractMatrix}}

    # add optional third argument for the case when the second set of weights is different.
    # Need to add a method for computing p_\bm too.

    # initialise stuff
    S = [([0],[0])] # needs to be a collection of pairs of indices
    f = [false] # S hasn't been processed yet.
    R = Tuple{Vector{Int64}, Vector{Int64}}[] # blank version of S
    μ₁ = getweights(Γ)
    G₁ = get_symmetry_group(Γ)
    M = length(Γ.IFS)
    A = zeros(1,1)
    B = zeros(1,1)
    L = zeros(1) # constant vector of log terms, only non-zero when s=0

    μ₁ == μ₂ ? fubuni_flag = true : fubuni_flag = false

    function scaler(ρ::Float64, m::Vector{<:Int64},m_::Vector{<:Int64},n::Vector{<:Int64},n_::Vector{<:Int64})
        # account for convention Γ₀:=Γ
        m  != [0] ? pₘ = prod(μ₁[m]) : pₘ = 1.0
        m_ != [0] ? pₘ_ = prod(μ₂[m_]) : pₘ_ = 1.0
        return ρ^(-s)*pₘ*pₘ_/prod(μ₁[n])/prod(μ₂[n_])
    end

    A_rows = 0
    A_cols = 0
    B_rows = 0
    B_cols = 0

    while sum(f)<length(f) # while f contains zeros
        for ∫∫_count = 1:length(S)
            if ~f[∫∫_count]
                a_row = zeros(length(S))
                a_row[∫∫_count] = 1.0
                b_row = zeros(length(R))
                ∫∫_indices = S[∫∫_count]
                for m=1:M, m_=1:M
                    mcat = vcat_(∫∫_indices[1],m)
                    mcat_ = vcat_(∫∫_indices[2],m_)
                    is_S_similar, ρ, similar_index = check_for_similar_integrals(Γ, S, mcat, mcat_, G₁, G₂, fubuni_flag)
                    is_ℓ_singular = check_for_ℓ_singular_integrals(Γ, mcat, mcat_)

                    # only need to check for R similarities if regular integral.
                    is_singular = is_S_similar || is_ℓ_singular
                    if !is_singular
                        is_R_similar, ρ, similar_index = check_for_similar_integrals(Γ, R, mcat, mcat_, G₁, G₂, fubuni_flag)
                    else
                        is_R_similar = false
                    end
                    # compute (3) from Dave's notes:
                    is_similar = is_S_similar || is_R_similar
                    if is_similar
                        is_S_similar ? similar_indices = S[similar_index] : similar_indices = R[similar_index]
                        # scale_adjust = 1/scaler(ρ, ∫∫_indices[1], ∫∫_indices[2], mcat, mcat_)
                        scale_adjust = 1/scaler(ρ, similar_indices[1], similar_indices[2], mcat, mcat_)
                    end

                    if is_ℓ_singular && !is_S_similar # new singularity type
                        push!(S,(mcat,mcat_)) # new type of singularitiy
                        push!(f,false)
                        push!(a_row, -1.0) # increase row by one
                        push!(L,0.0)
                        if A_rows > 0
                            A = hcat(A, zeros(A_rows))
                        end
                    elseif is_S_similar # singular, but seen similar
                        a_row[similar_index] -= scale_adjust
                        if s == 0
                             L[∫∫_count] += Γ.measure^2*log(1/ρ)*prod(μ₁[mcat])*prod(μ₂[mcat_]) # log constant adjustment
                        end
                    elseif is_R_similar # smooth, but seen similar
                        b_row[similar_index] += scale_adjust
                        if s == 0
                            L[∫∫_count] += Γ.measure^2*log(1/ρ)*prod(μ₁[mcat])*prod(μ₂[mcat_]) # log constant adjustment
                        end
                    else # smooth, nothing similar
                        push!(R,(mcat,mcat_))
                        push!(b_row, 1.0) # increase row by one
                        if B_rows > 0
                            B = hcat(B, zeros(B_rows))
                        end
                    end
                end
                # add new row to each matrix
                if A_rows == 0
                    A = reshape(a_row, 1, length(a_row))
                else
                    A = vcat(A, reshape(a_row, 1, length(a_row)))
                end
                if B_rows == 0
                    B = reshape(b_row, 1, length(b_row))
                else
                    B = vcat(B, reshape(b_row, 1, length(b_row)))
                end
                f[∫∫_count] = true
            end
            # update matrix sizes
            A_rows, A_cols = size(A)
            B_rows, B_cols = size(B)
        end
    end
    return A,B,S,R,L
end

"""
    s_energy(Γ::SelfSimilarFractal, s::Number, quad_rule::Function; μ₂::Vector{Float64} = getweights(Γ),
     G::Vector{AutomorphicMap}=TrivialGroup(Γ.spatial_dimension),
     G₁::Vector{AutomorphicMap}=TrivialGroup(Γ.spatial_dimension), G₂::Vector{AutomorphicMap}=TrivialGroup(Γ.spatial_dimension))


s is the value in |x-y|⁻ˢ, unless s==0, in which case log|x-y| is used.
μ₂ is an optional set of (probability) weights describing an invariant measure of the outer integral.
G₁ and G₂ are groups describing the symmetries of the inner and outer measures respectively.
If G is defined, both measures are assigned this symmetry.
Computes the s-energy of a fractal Γ, using the function quad_rule. This must be of the form:

    quad_rule = (e,j,f) -> I ≈ ∫ₑ∫ⱼ f(x,y) μ₁(x)μ₂(y)

where A and B are SelfSimilarFractal.
If quad_rule is replaced by some h::Number, the barycentre rule is used with meshwidth h.
"""
function s_energy(Γ::SelfSimilarFractal{V,M}, s::Number, ∫∫::Function;
                μ₂::Vector{Float64} = getweights(Γ),
                G₂::Vector{AutomorphicMap{V,M}} = TrivialGroup(Γ.spatial_dimension)) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}

    if getweights(Γ) == μ₂
        G₂=get_symmetry_group(Γ)
        Γ_μ₂ = Γ
    else
        Γ_μ₂ = changeweights(Γ,μ₂)
    end

    A,B,_,R,L = construct_singularity_matrix(Γ, s, μ₂=μ₂, G₂=G₂)
    
    r = zeros(length(R))
    for n=1:length(r)
        (m,m_) = R[n]
        # x,y,w = quad_rule(Γ[m],Γ_μ₂[m_])
        # r[n] = w'*Φₜ.(s,x,y)
        r[n] = ∫∫(Γ[m],Γ_μ₂[m_],(x,y)->Φₜ(s,x,y))
    end

    x = A\(B*r+L)

    return x[1]
end

# default to barycentre rule as follows: 
function s_energy(Γ::SelfSimilarFractal{V,M}, s::Number, h::Real; μ₂::Vector{Float64}=getweights(Γ),
     G₂::Vector{AutomorphicMap{V,M}} = TrivialGroup(Γ.spatial_dimension)) where {V<:Union{Real,AbstractVector}, M<:Union{Real,AbstractMatrix}}
    return  s_energy(Γ, s, (A::SelfSimilarFractal{V,M}, B::SelfSimilarFractal{V,M}, f::Function)->long_bary(A,B,f,h); μ₂ = μ₂, G₂=G₂)
end