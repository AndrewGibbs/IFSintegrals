# The below variation of vcat is needed for our convention that Γ₀:=Γ
vcat_(x,y) = vcat(x[x .!= 0],y[y .!= 0])
hcat_(x,y) = hcat(x[x .!= 0],y[y .!= 0])

function check_if_similar(Γ::SelfSimilarFractal,m,n,m_,n_, G, G_)
    # get shorthand for IFS
    S = Γ.IFS
    test_one = false
    test_two = false
    test_three = false
    is_similar = false # requires above three to be true
    ρ = 0.0 # initialise

    # define identity similarity, which is a workaround for index [0]
    s₀ = Similarity(1.0,zeros(Γ.spatial_dimension))

    m != [0] ? sₘ = sim_comp(S,m) : sₘ = s₀
    n != [0] ? sₙ = sim_comp(S,n) : sₙ = s₀
    m_ !=[0] ? sₘ_ = sim_comp(S,m_) : sₘ_ = s₀
    n_ !=[0] ? sₙ_ = sim_comp(S,n_) : sₙ_ = s₀

    # first test (6)
    if isapprox(sₘ.r/sₙ.r, sₘ_.r/sₙ_.r, atol=100*eps())
        ρ = sₘ.r/sₙ.r
        test_one = true
    end

    for T ∈ G, T_ ∈ G_
        # second test (7)
        if isapprox(sₘ.A*T.A*inv(sₙ.A) , sₘ_.A*T_.A*inv(sₙ_.A), atol=100*eps())
            test_two = true
        else
            test_two = false
        end
        
        # third test (8)
        if isapprox(sₘ.δ .- sₘ_.δ .- sₘ.r*sₘ.A*(T.A/sₙ.r*inv(sₙ.A)*sₙ.δ .- T.δ), .- sₘ_.r*sₘ_.A*(T_.A/sₙ_.r*inv(sₙ_.A)*sₙ_.δ .- T_.δ), atol=100*eps())
            #norm(sₘ.δ - sₘ_.δ - sₘ.r*sₘ.A*(T.A/sₙ.r*inv(sₙ.A)*sₙ.δ-T.δ) + sₘ_.r*sₘ_.A*(T_.A/sₙ_.r*inv(sₙ_.A)*sₙ_.δ - T_.δ)) ≈ 0
            test_three = true
        else
            # if test_two && test_one
            #     println(norm(sₘ.δ .- sₘ_.δ .- sₘ.r*sₘ.A*(T.A/sₙ.r*inv(sₙ.A)*sₙ.δ .- T.δ) .+ sₘ_.r*sₘ_.A*(T_.A/sₙ_.r*inv(sₙ_.A)*sₙ_.δ .- T_.δ)))
            # end
            test_three = false
        end

        # if tests are satisfies for some particular group element, can break loop early:
        if test_two && test_three
            break
        end
    end
    is_similar = test_one && test_two && test_three
    # if is_similar && abs(length(m)-length(n))==2 && m!=[0] && n!=[0]
    #     println(m,m_)
    #     println(n,n_)
    #     println("")
    # end
    return is_similar, ρ
end

function check_for_similar_integrals(Γ, X,mcat, mcat_, G₁, G₂, fubini_flag) # here X is some set of indices... leave abstract for now
    is_X_similar = false
    similar_index = nothing
    proportionality_const = 0.0
    # should compactly write the following as a function, it's almost repeated
    for ∫∫_index = 1:length(X)
        ∫∫_indices_higher_level = X[∫∫_index]
        # check for opposite index, i.e. (m,m_) ∼ (m_,m), assuming μ₁ = μ₂
        # if fubini_flag && X[∫∫_index][1] == mcat_ && X[∫∫_index][2] == mcat
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

function check_for_ℓ_singular_integrals(Γ::SelfSimilarFractal, mcat, mcat_)
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

function construct_singularity_matrix(Γ::SelfSimilarFractal, s::Number; μ₂::Vector{Float64} = getweights(Γ), G₁=TrivialGroup(Γ.spatial_dimension), G₂=TrivialGroup(Γ.spatial_dimension))

    # add optional third argument for the case when the second set of weights is different.
    # Need to add a method for computing p_\bm too.

    # initialise stuff
    S = [([0],[0])] # needs to be a collection of pairs of indices
    f = [false] # S hasn't been processed yet.
    R = Tuple{Vector{Int64}, Vector{Int64}}[] # blank version of S
    μ₁ = Γ.weights
    # B = Matrix{Float64}[]
    M = length(Γ.IFS)
    A = zeros(1,1)
    B = zeros(1,1)

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
                        # print((mcat, mcat_))
                        # print(",")
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
                        if A_rows > 0
                            A = hcat(A, zeros(A_rows))
                        end
                    elseif is_S_similar # singular, but seen similar
                        a_row[similar_index] -= scale_adjust
                    elseif is_R_similar # smooth, but seen similar
                        b_row[similar_index] += scale_adjust
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
    return A,B,S,R
end

function s_energy(Γ::SelfSimilarFractal, s::Number, quad_rule::Function; μ₂::Vector{Float64} = getweights(Γ), G=nothing, G₁=TrivialGroup(Γ.spatial_dimension), G₂=TrivialGroup(Γ.spatial_dimension))

    if G!==nothing
        G₁ = G
        G₂ = G
    end

    A,B,_,R = construct_singularity_matrix(Γ, s, μ₂=μ₂, G₁=G₁, G₂=G₂)

    μ₁ = Γ.weights
    if μ₁ == μ₂
        Γ_μ₂ = Γ
    else
        Γ_μ₂ = changeweights(Γ,μ₂)
    end
    
    r = zeros(length(R))
    # num_pts = 0.0
    for n=1:length(r)
        (m,m_) = R[n]
        x,y,w = quad_rule(Γ[m],Γ_μ₂[m_])
        # num_pts += length(w)
        r[n] = w'*Φₜ.(s,x,y)
    end
    # println(num_pts)

    x = A\(B*r)

    return x[1]
end

# default to barycentre rule as follows:
s_energy(Γ::SelfSimilarFractal, s::Number, h::Real; μ₂::Vector{Float64}=getweights(Γ), G=nothing, G₁=TrivialGroup(Γ.spatial_dimension), G₂=TrivialGroup(Γ.spatial_dimension)) = s_energy(Γ, s, (A::SelfSimilarFractal, B::SelfSimilarFractal)->barycentre_rule(A,B,h); μ₂ = μ₂, G=G, G₁=G₁, G₂=G₂)
# s_energy(Γ::SelfSimilarFractal, s::Number, h::Real; μ₂::Vector{Float64} = getweights(Γ), G=nothing) = s_energy(Γ,s,(A::SelfSimilarFractal, B::SelfSimilarFractal)->barycentre_rule(A,B,h); μ₂=μ₂, G=nothing)