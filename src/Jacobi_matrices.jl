
function get_modified_coeffs_from_coeffs(coeffs_one_below::Vector{Float64},coeffs_two_below::Vector{Float64},A::Vector{Float64},r::Vector{Float64},sₘ::Similarity)
    # rule of thumb: square bracket indices should always be incremented by one compared with what's written above,
     # to account for the zero index

     # Conflicting notation between Mantica '96 and myself, so this step should switch everything to his notation:
     δᵢ = sₘ.r
     βᵢ = sₘ.δ
     # that's better.
     
     n = length(coeffs_one_below) # one more because of the zeroth entry, one less because we're incrementing
     modified_coeffs = zeros(Float64,n+1)
     
     function γ(ℓ::Integer,j::Integer)
         if ℓ==j
             return A[ℓ+1]*coeffs_one_below[j]
         elseif ℓ==j-1
             return r[ℓ+1]*coeffs_one_below[j]
         elseif j==ℓ-1
             return r[j+1]*coeffs_one_below[j]
         else
             return 0.0
         end
     end

     
# Currently, you're mixing the modified coefficients for $\tilde{p}$ with the (normalised) coefficients for $p$.
# Need to write a function which maps between these, given $r$.

         
     function get_trickier_sum(j)
         trickier_sum = 0.0
         for ℓ=max(0,j-1):min(n-1,j+1)
             trickier_sum += γ(ℓ,j)
         end
         return trickier_sum
     end
         
     for j = 0:(n-2)
         modified_coeffs[j+1] = (βᵢ-A[n])*coeffs_one_below[j+1] +  δᵢ*get_trickier_sum(j+1) - r[n]*coeffs_two_below[j+1]
     end
     
     # now deal with j = n-1
     modified_coeffs[n] = (βᵢ-A[n])*coeffs_one_below[n] +  δᵢ*get_trickier_sum(n)
     
     # now deal with j = n, which is not in orthanormal function
     modified_coeffs[n+1] =  δᵢ*coeffs_one_below[n]
     
     return modified_coeffs
 end

 function get_rₙ(modified_coeffs::Vector{Vector{Float64}}, coeffs_one_below::Vector{Vector{Float64}}, A::Vector{Float64}, Γ::SelfSimilarFractal)
    n = length(modified_coeffs::Vector{Vector{Float64}})-1
     M = length(Γ.IFS)
     B = zeros(Float64,M)
     C = zeros(Float64,M)
     D = zeros(Float64,M)
     D_nominator = zeros(Float64,M)
     
     for i=1:M
         for ℓ=0:(n-1)
            B[i] += (Γ.IFS[i].δ + Γ.IFS[i].r*A[ℓ+1])*modified_coeffs[i][ℓ+1]*coeffs_one_below[i][ℓ+1]
         end
         for ℓ=0:(n-2)
            C[i] +=  Γ.IFS[i].r*r[ℓ+2]*(modified_coeffs[i][ℓ+1]*coeffs_one_below[i][ℓ+2] + modified_coeffs[i][ℓ+2]*coeffs_one_below[i][ℓ+1])
         end
         D_nominator[i] = Γ.IFS[i].r*modified_coeffs[i][n+1]*coeffs_one_below[i][n]
     end
     
     rₙ² = (Γ.weights'*(B+C)) / (1-(Γ.weights'*D_nominator))
     
     return sqrt(rₙ²)
 end

 function get_Aₙ(coeffs_this_level::Vector{Vector{Float64}}, A::Vector{Float64}, r::Vector{Float64}, Γ::SelfSimilarFractal)
    M = length(Γ.IFS)
    n = length(coeffs_this_level)-1
    big_sum = 0.0
    denominator = 1.0
    for i=1:M
        big_sum += coeffs_this_level[i][n+1]^2*Γ.IFS[i].δ
        for m=0:(n-1)
            big_sum += coeffs_this_level[i][m+1]^2*(Γ.IFS[i].δ+Γ.IFS[i].r*A[m+1]) + coeffs_this_level[i][m+1]*coeffs_this_level[i][m+2]*Γ.IFS[i].r*(r[m+1]+r[m+2])
        end
        big_sum *= Γ.weights[i]
        denominator += Γ.weights[i]*coeffs_this_level[i][n+1]^2*Γ.IFS[i].r
    end
    return big_sum/(1-denominator)
end

function get_Jacobi_matrix(Γ::Attractor,N::Int64)
    # initialisation
   A = [Γ.measure*Γ.barycentre] # checks out, given def'n of barycentre, should be ∫_Γ x dμ(x)
    r = [0.0]
    M = length(Γ.IFS)
    coeffs_one_below = [[1.0] for _=1:M]
    coeffs_two_below = [[0.0] for _=1:M]
    
    # iteration
    for n=1:N #we've done n=0 above
        # initialise this level
        coeffs_this_level = [zeros(n+1) for _=1:M]
        
        # step one
        modified_coeffs = [zeros(n+1) for _=1:M]
        for m=1:M
            modified_coeffs[m] = get_modified_coeffs_from_coeffs(coeffs_one_below[m], coeffs_two_below[m], A, r, Γ.IFS[m])
        end
        
        # step two
        rₙ = get_rₙ(modified_coeffs, coeffs_one_below, A, Γ)
        push!(r,rₙ)
        
        # step three
        for m=1:M
            coeffs_this_level[m][1:n] = modified_coeffs[m][1:n]/r[n+1]
            coeffs_this_level[m][n+1] = modified_coeffs[m][n+1]
        end
        
        # step four
        Aₙ = get_Aₙ(coeffs_this_level, A, r, Γ)
        push!(A,Aₙ)
        
        coeffs_two_below = coeffs_one_below
        coeffs_one_below = coeffs_this_level
        
    end
    
    # now make the Jacobi matrix
    J = zeros(Float64,N+1,N+1)
    for i=0:(N-1)
       J[i+1,i+1] = A[i+1]
        J[i+1,i+2] = r[i+2]
        J[i+2,i+1] = r[i+2]
    end
    J[N+1,N+1] = A[N+1]
    
    return J
end