# The below variation of vcat is needed for our convention that Γ₀:=Γ
vcat_(x,y) = vcat(x[x .!= 0],y[y .!= 0])
vcat_(x,y) = hcat(x[x .!= 0],y[y .!= 0])

function check_if_similar(Γ::SelfSimilarFractal,m,n,m_,n_)
    # get shorthand for IFS
    S = Γ.IFS

    condition_two_holds = true
    ρ̃  = 0.0 # initialise

    sₘ = sim_comp(S,m)
    sₙ = sim_comp(S,n)
    sₘ_ = sim_comp(S,m_)
    sₙ_ = sim_comp(S,n_)

    # first test (6)
    if sₘ.r/sₙ.r ≈ sₘ_.r/sₙ_.r
        ρ̃ = sₘ.r/sₙ.r
    else
        condition_two_holds = false
    end
    
    # second test (7)
    if sₘ.A*inv(sₙ.A) ≈ sₘ_.A*inv(sₙ_.A)
        nothing
    else
        condition_two_holds = false
    end
    
    # third test (8)
    if norm(sₘ.δ - sₘ_.δ - sₘ.r/sₙ.r*sₘ.A*inv(sₙ.A)*sₙ.δ + sₘ_.r/sₙ_.r*sₘ_.A*inv(sₙ_.A)*sₙ_.δ) ≈ 0
        nothing
    else
        condition_two_holds = false
    end
        
    return [condition_two_holds, ρ̃ ]
end

function check_for_similar_integrals(X,mcat,mcat_) # here X is some set of indices... leave abstract for now
    is_X_similar == false
    similar_index = nothing
    # should compactly write the following as a function, it's almost repeated
    for ∫∫_index = 1:length(X)
        ∫∫_indices_higher_level = X[∫∫_index]
        this_is_X_similar, ρ̃  = check_for_condition_two(Γ,∫∫_indices_higher_level[1],mcat,_indices_higher_level[2],mcat_)
        if this_is_X_similar
            is_X_similar = true
            proportionality_const = ρ̃ 
            similar_index = ∫∫_index
            break
        end
    end
    return is_X_similar, proportionality_const, similar_index
end

function check_for_ℓ_singular_integrals(Γ_singularities::Matrix{Bool}, mcat, mcat_; ℓ_depth==2)
    is_singular = false

    if ℓ_depth != 2
        error("Haven't coded this yet")
    end
    
    if length(mcat)==1 && length(mcat_)==1
        Γ_singularities[mcat[1],mcat_[1]] ? is_singular = true : nothing
    end

    return is_singular
end

function construct_singularity_matrix(Γ::SelfSimilarFractal, Γ_singularities::Matrix{Bool}, s::Float64; μ₂::Vector{Float64} = [0.0])

    # add optional third argument for the case when the second set of weights is different.
    # Need to add a method for computing p_\bm too.

    if isa(Γ,Attractor)
        μ₁ = Γ.weights
    else
        μ₁ = Γ.attractor.weights
    end
    if μ₂ == [0.0]
        μ₂ = μ₁
    end
    for μ = [μ₁,μ₂]
        length(μ)!=M ? error("μ must be a vector containing the same length as the IFS") : nothing
        sum(μ) ≈ 1 ? nothing : error("μ must sum to one")
    end

    function scaler(ρ::Float64, m::Vector{<:Int64},m_::Vector{<:Int64},n::Vector{<:Int64},n_::Vector{<:Int64})
        # account for convention Γ₀:=Γ
        m == [0] ? pₘ = prod(μ₁[m]) : pₘ = 1.0
        m_ == [0] ? pₘ_ = prod(μ₂[m_]) : pₘ_ = 1.0
        return ρ^(-s)*pₘ*pₘ_/prod(μ₁[n])/prod(μ₂[n_])
    end

    # initialise stuff
    S = [([0],[0])] # needs to be a collection of pairs of indices
    f = [false] # S hasn't been processed yet. 
    R = Vector{Tuple{Vector{Int64}, Vector{Int64}}}[] # blank version of S
    # B = Matrix{Float64}[]
    M = length(Γ.IFS)

    A_rows = 0
    A_cols = 0
    B_rows = 0
    B_cols = 0

    while sum(f)<length(f) # while f contains zeros
        for ∫∫_count = 1:length(S)
            if ~f[∫∫_count]
                a_row = [0.0]
                b_row = [0.0]
                ∫∫_indices = S[∫∫_count]
                for m=1:M, m_=1:M
                    mcat = vcat_(∫∫_indices[1],m)
                    mcat_ = vcat_(∫∫_indices[2],m_)

                    is_S_similar, ρ, similar_index = check_for_similar_integrals(S, mcat, mcat_)
                    is_ℓ_singular = check_for_ℓ_singular_integrals(Γ_singularities, mcat, mcat_)

                    # only need to check for R similarities if regular integral.
                    is_singular = is_S_similar || is_ℓ_singular
                    if !is_singular
                        is_R_similar, ρ, similar_index = check_for_similar_integrals(R, mcat, mcat_)
                    else
                        is_R_similar = false
                    end
                    # compute (3) from Dave's notes:
                    scale_adjust = scaler(ρ, ∫∫_indices[1], ∫∫_indices[2], mcat, mcat_)

                    if is_ℓ_singular && !is_S_similar # new singularity type
                        push!(S,(mcat,mcat_)) # new type of singularitiy
                        push!(f,false)
                        push!(a_row,0.0) # increase row by one
                        if A_rows > 0
                            A = hcat(A, zeros(A_rows))
                        end
                    elseif is_S_similar # singular, but seen similar
                        a_row[similar_index] *= scale_adjust
                    elseif is_R_similar # smooth, but seen similar
                        b_row[similar_index] *= scale_adjust
                    else # smooth, nothing similar
                        push!(R,(mcat,mcat_))
                        push!(b_row,0.0) # increase row by one
                        if B_rows > 0
                            B = hcat(B, zeros(B_rows))
                        end
                    end
                end
                # add new row to each matrix
                if A_rows == 0
                    A = reshape(a_row, 1, A_cols)
                else
                    A = vcat(A, a_row)
                end
                if B_rows == 0
                    B = reshape(b_row, 1, B_cols)
                else
                    B = vcat(B, b_row)
                end
                f[int_count] = true
            end
        end
    end
    return A,B,S,R
end