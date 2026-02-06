# Last revision: 2025-12-10 19:12

function synth_cnv(fragment_coords, fragment_status, threshold_proportion)
    # CNV synthesis function: merge CNV fragments into CNV regions
    # Parameters:
    # - fragment_coords: n×4 data frame
    #     Column 1: Chromosome identifier (string)
    #     Column 2: Start position of current scaffold fragment
    #     Column 3: End position of current scaffold fragment
    #     Column 4: Length of current scaffold fragment
    # - fragment_status: numeric vector of length n with values {0, 1, -1}
    #     0: current scaffold fragment is negative
    #     1: current scaffold fragment is positive
    #     -1: current scaffold fragment has undetermined status

    # Return value:
    # - m×2 matrix
    #     Column 1: Row index of CNV region start fragment in input matrix
    #     Column 2: Row index of CNV region end fragment in input matrix
    # - Returns -1 if an error occurs


    # Step 1: Check if parameter 1 is sorted by chromosome and start position
    # Check if chromosomes are ordered
    chromosomes = fragment_coords[:, 1]
    start_positions = fragment_coords[:, 2]
    
    # Check sorting by chromosome and start position
    for i in 2:size(fragment_coords, 1)
        # If chromosome number decreases, not sorted
        if chromosomes[i] < chromosomes[i-1]
            println("Error: CNV fragment coordinate matrix is not sorted by chromosome and start position")
            return -1
        # If on same chromosome but start position decreases, not sorted
        elseif chromosomes[i] == chromosomes[i-1] && start_positions[i] < start_positions[i-1]
            println("Error: CNV fragment coordinate matrix is not sorted by chromosome and start position")
            return -1
        end
    end
    
    # Step 2: Select row indices from parameter 1 (scaffold fragment positions) where parameter 2 (copy difference status) equals 1 to form CNV table
    # Get row indices of positive fragments (value = 1)
    positive_indices = findall(fragment_status .== 1)
    
    # If no positive fragments, return empty matrix
    if isempty(positive_indices)
        println("Error: No positive fragments found")
        return -1
    end

    # Step 3: Merge CNV regions based on positive fragment content after proposed merging
    # Initialize CNV table with both columns containing positive fragment row indices
    #= Use a k×2 matrix to simulate a linked list with deletion capability for recording fragment merging.
    Column 1: fragment start position, Column 2: fragment end position.
    When merging fragments i and i+1 (deletion operation), set row i column 2 to -1, and row i+1 column 1 to -1.
    For the entire matrix, the j-th non-(-1) element in column 1 always pairs with the j-th non-(-1) element in column 2 logically as (start_pos, end_pos), recording a fragment.
    When outputting final result, reorganize by writing all non-(-1) elements to corresponding columns to form m×2 output matrix.
    =#
    cnv_table = hcat(positive_indices, positive_indices)
    
    k = size(cnv_table, 1)
    if k <= 1
        # Only one CNV fragment, no merging needed
        return cnv_table
    end    

    # Evaluate positive fragment content after proposed merging to guide merging
    # In each evaluation round, merge all that meet requirements
    # Perform multiple evaluation rounds until no more merging is possible (no fragments meet requirements)

    merged = true
    while merged # Repeatedly evaluate until no more merging is possible
        merged = false
        # Start one round of evaluation
        for i in 2:size(cnv_table, 1) 
            # Find fragment indices represented by these two rows in original CNV table
            last_idx_fragment_coords = cnv_table[i-1, 1]
            current_idx_fragment_coords = cnv_table[i, 2]            

            # Calculate positive fragment content
            temp_lengths = fragment_coords[last_idx_fragment_coords:current_idx_fragment_coords, 4]
            possible_proportion = sum(temp_lengths .* (fragment_status[last_idx_fragment_coords:current_idx_fragment_coords] .== 1)) / 
                                  sum(temp_lengths)           

            if (possible_proportion >= threshold_proportion) & 
                (fragment_coords[last_idx_fragment_coords, 1] == fragment_coords[current_idx_fragment_coords, 1])
                # Meets merging requirements, perform merge
                cnv_table[i, 1] = cnv_table[i - 1, 1] # Merge
                cnv_table[i - 1, 1] = -1 # Delete
                cnv_table[i - 1, 2] = -1 # Delete
                merged = true # If merged, need next evaluation round
            end            
        end
        # After each evaluation and merging round, update cnv_table
        temp_cnvtable = hcat(cnv_table[findall(cnv_table[:, 1] .> 0), 1], cnv_table[findall(cnv_table[:, 2] .> 0), 2])
        cnv_table = temp_cnvtable
    end     
    
    return cnv_table
end


function upgrade_state_vector(prob_vector, 
                             state_vector, 
                             fdr_threshold)
# Parameters
# - `prob_vector::Vector{Float64}`: Probability vector, each element is a probability value in [0,1]
# - `state_vector::Vector{Int}`: State vector, values are 0 (negative) or 1 (positive)
# - `fdr_threshold::Float64`: FDR control threshold, e.g., 0.05
# Return
# - `Vector{Int}`: Upgraded state vector
# Workflow
# 1. Identify indices of negative positions adjacent to positive states
# 2. Extract probability values at corresponding positions
# 3. Perform FDR correction using Benjamini-Hochberg method
# 4. Update state vector based on significance results
    
    # Parameter validation
    @assert length(prob_vector) == length(state_vector) # "Probability and state vectors must have same length"
    @assert all(0 .<= prob_vector .<= 1) # "Probability values must be in [0,1] range"
    @assert all(state_vector .== 0 .|| state_vector .== 1) # "State vector values must be 0 or 1"
    @assert 0 <= fdr_threshold <= 1 # "FDR threshold must be in [0,1] range"

    
    
    n = length(state_vector)
    adjacent_zero_indices = Int[]
    
    # Step 1: Find negative positions adjacent to positive states
    for i in 1:n
        if state_vector[i] == 1
            # Check left neighbor
            if i > 1 && state_vector[i-1] == 0
                push!(adjacent_zero_indices, i-1)
            end
            # Check right neighbor
            if i < n && state_vector[i+1] == 0
                push!(adjacent_zero_indices, i+1)
            end
        end
    end

    

    # Remove duplicate indices
    adjacent_zero_indices = unique(adjacent_zero_indices)
    #println("Length of adjacent negative position indices: ", length(adjacent_zero_indices))
   
    
    # Step 2: Extract probability values at adjacent negative positions
    if isempty(adjacent_zero_indices)
        return state_vector  # Return directly when no adjacent negative states
    end
    
    adjacent_probs = prob_vector[adjacent_zero_indices]    

    # Step 3: FDR correction and significance determination
    adjusted_pvalues = adjust(adjacent_probs, BenjaminiHochberg())
    
    significant_indices = findall(adjusted_pvalues .< fdr_threshold)
    #println("Length of significant adjacent negative indices: ", length(significant_indices)) 
    #println("Sum of adjacent significant negative probabilities: ",  -log10.(adjacent_probs'))
    #println("Difference: ", sum((state_vector[adjacent_zero_indices]  .== 0)))

    if isempty(significant_indices)
        return state_vector  # Return directly when no adjacent negative states meet threshold
    end
    
    # Update state vector
    for i in 1:length(significant_indices)
        original_idx = adjacent_zero_indices[significant_indices[i]]
        state_vector[original_idx] = 1
    end
    
    # Step 4: Return upgraded state vector
    return state_vector
end


#***********************Enhanced copy number variation test****************************
function prob_copy_differenceM(frag_depths, num_copies, givenParameters, windowSize)
    # frag_depth: n×2 matrix, each column represents sequencing depth of short fragments for sample 1 and sample 2 respectively
    # num_copies: copy numbers for sample 1 and sample 2 of this fragment
    # givenParameters: estimated parameters, return values of estimate_gmm() for sample 1 and sample 2
    # windowSize: fragment length used when estimating givenParameters
    # Returns: probability that no copy number variation exists between the two samples

    # Dependency: Distributions.jl

    # Number of components
    len_frag = size(frag_depths)[1]
    num_copies = Int64.(num_copies)
    #pesudo_len = len_frag/windowSize
    pesudo_len = sqrt(len_frag/windowSize)
    #pesudo_len = 1

    meanwithin =  mean(frag_depths, dims = 1) # 1×2 matrix, columns are different samples
    meantotal = [givenParameters[1].means[num_copies]; givenParameters[2].means[num_copies]] # 2×2 matrix, rows are different fragments, columns are different samples
    varwithin =  [var(frag_depths[:,1]) var(frag_depths[:,2])] ./ windowSize # 1×2 matrix, columns are different samples
    vartotal = [givenParameters[1].stds[num_copies] .^ 2; givenParameters[2].stds[num_copies] .^ 2] # 2×2 matrix, rows are different fragments, columns are different samples

    # varbetween_s1 = max.(10^-16, reshape(vartotal[1,:],1,:) .- varwithin) # Between-group variance, must be > 0
    # varbetween_s2 = max.(10^-16, reshape(vartotal[2,:],1,:) .- varwithin)

    varbetween_s1 = max.(10^-16, reshape(vartotal[1,:],1,:) .- varwithin) # 1×2 matrix, between-group variance for fragment 1 in different groups, must be > 0
    varbetween_s2 = max.(10^-16, reshape(vartotal[2,:],1,:) .- varwithin) # 1×2 matrix, between-group variance for fragment 2 in different groups, must be > 0

    chi_dist_s1 = (meanwithin .- meantotal[1, :]) .^ 2 ./ (varwithin ./ pesudo_len + varbetween_s1) # 1×2 matrix, chi-square values for fragment 1 in different sample groups
    chi_dist_s2 = (meanwithin .- meantotal[2, :]) .^ 2 ./ (varwithin ./ pesudo_len + varbetween_s2) # 1×2 matrix, chi-square values for fragment 2 in different sample groups


    chi_dist = Chisq(2)       # Define central chi-square distribution with 2 degrees of freedom
    p_value_g1 = 1 - cdf(chi_dist, sum(chi_dist_s1)) 
    p_value_g2 = 1 - cdf(chi_dist, sum(chi_dist_s2)) 

    output = (p_value_g1 * givenParameters[1].weights[num_copies[1]] + p_value_g2 * givenParameters[2].weights[num_copies[2]]) /
             (givenParameters[1].weights[num_copies[1]] + givenParameters[2].weights[num_copies[2]])

    return output
end

#*****************Functions below are for generating Markov chains******************************
function simulate_markov_chain(
    M::Matrix{Float64},  # Transition probability matrix (n×n)
    C::Int,              # Number of time steps
    initial_state::Int   # Initial state (1 to n)
)::Vector{Int}
    # This function generates a vector of length C based on transition probability matrix M
    # Elements of output vector chain represent the state at each time point
    # Output states are between 1 and n; consider whether to adjust to 0 to n-1 for CNV simulation
    # Dependency: "StatsBase"
    
    # Parameter validation
    @assert size(M,1) == size(M,2) "Transition matrix must be square"
    @assert all(0 .<= M .<= 1) "Probability values must be in [0,1] range"
    @assert all(sum(M, dims=2) .≈ 1) "Each row probability sum must equal 1"
    @assert 1 <= initial_state <= size(M,1) "Initial state out of bounds"
    
    states = Vector{Int}(undef, C)
    states[1] = initial_state 
    
    # State transition simulation
    for t in 2:C
        states[t] = sample(1:size(M,2), Weights(M[states[t-1], :]))
    end

    return states
end

#*****************Functions below are for generating confounding factors like GC content, mapping errors, etc.******************************
function generate_multinomial_samples(
    p::Vector{Float64},  # Probability vector (probability of each event)
    n_trials::Int,       # Number of trials
    sample_size::Int     # Number of samples to generate
)::Vector{Vector{Int}}
    # Dependency: "Distributions"

    # Parameter validation
    @assert all(0 .<= p .<= 1) "Probability values must be in [0,1] range"
    @assert sum(p) ≈ 1 "Probability vector sum must equal 1"
    @assert n_trials > 0 "Number of trials must be positive integer"
    
    # Create multinomial distribution
    d = Multinomial(n_trials, p)
    
    # Generate samples (each sample is a vector of length equal to length(p))
    return [rand(d) for _ in 1:sample_size]
end

#****************************************************************************************
function estimate_gmm(data, n_components; means=0, max_iter=400, tol=1e-6)
    # Estimate parameters of Gaussian mixture model
    # This version eliminates 'division by zero' errors caused by insufficient numerical precision

    # Initialize parameters
    n = length(data)
    weights = fill(1/n_components, n_components)
    if length(means) != n_components
          means = quantile(data, range(0, 1, length=n_components+2)[2:end-1])
    end
    #means = quantile(data, range(0, 1, length=n_components+2)[2:end-1])
    stds = fill(std(data)/n_components, n_components)
    
    # EM algorithm
    for iter in 1:max_iter        
        # E-step: Calculate posterior probabilities
        responsibilities = zeros(n, n_components)
        for k in 1:n_components            
            responsibilities[:,k] = weights[k] * pdf.(Normal(means[k], stds[k]), data)
            temp_posi = findall(responsibilities[:,k] .< 10^-300) 
            responsibilities[temp_posi,k] .= 10^-300  # Limit minimum value to avoid 'division by zero' error
        end
        responsibilities ./= sum(responsibilities, dims=2) 
        # M-step: Update parameters
        new_means = zeros(n_components)
        new_stds = zeros(n_components)
        for k in 1:n_components
            resp_k = responsibilities[:,k]            
            new_means[k] = sum(resp_k .* data) / sum(resp_k)
            new_stds[k] = sqrt(sum(resp_k .* (data .- new_means[k]).^2) / sum(resp_k))
            weights[k] = mean(resp_k)
        end
        
        # Check convergence
        if maximum(abs.(new_means - means)) < tol && 
           maximum(abs.(new_stds - stds)) < tol
            break
        end        

        if iter == max_iter
            max_dev_mean = maximum(abs.(new_means - means))
            max_dev_std = maximum(abs.(new_stds - stds))
            println("Maximum iteration reached at iter = $iter")
			println("$max_dev_mean")
			println("$max_dev_std")
			
        end
        
        means, stds = new_means, new_stds
    end
    
    return (means=means, stds=stds, weights=weights)
end

#************************Calculate probability that left and right flanks have same copy number*************************************
function rate_coincident(obs_data, givenParameters)
    # Calculate probability matrix for copy number of two observed values
    # obs_data: vector of length 2, recording sequencing depth observations on left and right sides
    # givenParameters: tuple from function estimate_gmm
    # Dependency: Distributions.jl
    
    # Number of components
   #num_components = length(givenParameters.means)
   normal_obs_data = (obs_data' .-  givenParameters.means) ./ givenParameters.stds

   d = Normal(0, 1)       # Define standard normal distribution (μ=0, σ=1)
   #cdf_value = cdf(d, 1.96) 

   pdf_value = pdf(d, normal_obs_data) 

   pdr_value = (pdf_value .* givenParameters.weights) ./(sum(pdf_value .* givenParameters.weights, dims = 1))

   coincident_rate = sum(pdr_value[:, 1] .* pdr_value[:, 2])

   return coincident_rate
end

#***********************Count run lengths of non-zero values in a vector***************************************
function run_length_nonzero(v)
    # Initialize result matrix
    result = Matrix{Int}(undef, 0, 3)
    n = length(v)
    n == 0 && return result
    
    # Initialize run length variables
    current_val = v[1]
    start_idx = 1
    
    for i in 2:n
        if v[i] != current_val
            # When value changes, record previous run (if non-zero)
            if current_val != 0
                result = vcat(result, [current_val start_idx i-1])
            end
            current_val = v[i]
            start_idx = i
        end
    end
    
    # Handle the last run
    if current_val != 0
        result = vcat(result, [current_val start_idx n])
    end
    
    return result
end

#***********************Determine copy number of aggregated fragments***************************************
function identify_copy_number(frag_depth, givenParameters, threshold)
    # frag_depth: vector, copy number of each fragment
    # givenParameters: estimated parameters, return value of estimate_gmm()
    # threshold: threshold value
    # Returns: matrix where first column is copy number (-1 indicates undetermined), subsequent columns are probabilities for copy numbers 1, 2, 3, ...

    # Dependency: Distributions.jl

    # Number of components
    num_components = length(givenParameters.means)
    normal_obs_data = (frag_depth .-  givenParameters.means) ./ givenParameters.stds

    temp = zeros(num_components)

    d = Normal(0, 1)       # Define standard normal distribution (μ=0, σ=1)
    pdf_value = pdf(d, normal_obs_data) 
    pdr_value = (pdf_value .* givenParameters.weights) ./(sum(pdf_value .* givenParameters.weights, dims = 1))
    
    tempBayesFactors = pdr_value ./ (1 .- pdr_value)
    (max_val, max_idx) = findmax(tempBayesFactors)
    if max_val > threshold
        copyNum = max_idx
    else
        copyNum = -1
    end
    
    return [copyNum pdr_value']
end

#***********************Copy number variation test****************************
function prob_copy_difference(frag_depth, num_copies, givenParameters)
    # frag_depth: vector of length 2, sequencing depth of this fragment for sample 1 and sample 2
    # num_copies: copy numbers for sample 1 and sample 2 of this fragment
    # givenParameters: estimated parameters, return values of estimate_gmm() for sample 1 and sample 2
    # Returns: probability that no copy number variation exists between the two samples

    # Dependency: Distributions.jl

    # Number of components
    num_copies = Int64.(num_copies)
    normal_frag_depths_g1 = (frag_depth .-  givenParameters[1].means[num_copies[1]]) ./ givenParameters[1].stds[num_copies[1]]
    normal_frag_depths_g2 = (frag_depth .-  givenParameters[2].means[num_copies[2]]) ./ givenParameters[2].stds[num_copies[2]]

    chi_dist = Chisq(2)       # Define central chi-square distribution with 2 degrees of freedom
    p_value_g1 = 1 - cdf(chi_dist, sum(normal_frag_depths_g1.^2)) 
    p_value_g2 = 1 - cdf(chi_dist, sum(normal_frag_depths_g2.^2)) 

    output = (p_value_g1 * givenParameters[1].weights[num_copies[1]] + p_value_g2 * givenParameters[2].weights[num_copies[2]]) /
             (givenParameters[1].weights[num_copies[1]] + givenParameters[2].weights[num_copies[2]])

    return output
end

