
using Pkg

# Define list of required packages
required_packages = [
    "CSV", "DataFrames", "Plots", 
    "StatsBase", "Distributions", "DelimitedFiles", 
    "MultipleTesting" 
    ]

# Check and load each package
for pkg in required_packages
    try
        eval(Meta.parse("using $pkg"))
        println("$pkg loaded")
    catch
        println("$pkg not installed, installing...")
        Pkg.add(pkg)
        eval(Meta.parse("using $pkg"))
        println("$pkg installed and loaded successfully")
    end
end


"""
    process_genomic_data(temp_path::String, file_names::Vector{String})::DataFrame

Encapsulate the genomic data processing workflow into a single function.

This function reads multiple genomic data files, filters chromosomes, and combines 
the data into a single DataFrame with positional information and coverage data.

# Arguments
- `temp_path::String`: The path to the directory containing the data files
- `file_names::Vector{String}`: A vector of file names to process

# Returns
- `DataFrame`: A combined DataFrame containing:
    - Positional information (Chromosome, BinLowEdge, BinUpEdge)
    - Coverage data from multiple samples

# Examples
```julia
temp_path = "C:/Users/yghe_/Documents/Work/Julia/DuanAnQi/293细胞每1KbRD值/代码重构/"
file_names = ["E9-E10_his_rd_p_1000_allchr_new.csv", ...]
combineData = process_genomic_data(temp_path, file_names)
```

# Exceptions
- Throws ArgumentError if temp_path is empty
- Throws ArgumentError if file_names is empty
"""
function process_genomic_data(temp_path::String, file_names::Vector{String})::DataFrame
    # Input validation
    if isempty(temp_path)
        throw(ArgumentError("temp_path must be a non-empty string"))
    end
    
    if isempty(file_names)
        throw(ArgumentError("file_names must be a non-empty vector"))
    end

    # Read all data files
    dfs = Vector{DataFrame}(undef, length(file_names))
    for i in 1:length(file_names)
        full_path = joinpath(temp_path, file_names[i])
        dfs[i] = CSV.read(full_path, DataFrame)
    end

    # Define chromosome strings for filtering
    chr_strings = [
        "chr1", "chr2", "chr3", "chr4", "chr5", 
        "chr6", "chr7", "chr8", "chr9", "chr10", 
        "chr11", "chr12", "chr13", "chr14", "chr15", 
        "chr16", "chr17", "chr18", "chr19", "chr20", 
        "chr21", "chr22" 
    ]

    # Create flag for rows with defined chromosomes
    n_rows = size(dfs[1], 1)
    chr_flag = zeros(n_rows, 1)

    for i in 1:n_rows
        chr_flag[i] = (sum(dfs[1].Chromosome[i] .== chr_strings) > 0)
    end

    # Find indices with fully defined positions
    defined_indices = findall(chr_flag[:] .== 1)

    # Create combined DataFrame with positional information
    combineData = DataFrame(
        Chromosome = dfs[1][defined_indices, :Chromosome],
        BinLowEdge = dfs[1][defined_indices, :BinLowEdge],
        BinUpEdge = dfs[1][defined_indices, :BinUpEdge],
    )

    # Add data from other files
    for i in 1:length(dfs)
        # Generate column name from file name
        column_prefix = split(file_names[i], "_")[1]
        column_name = replace(column_prefix, "-" => "")
        column_symbol = Symbol("df" * column_name)
        
        combineData[!, column_symbol] = dfs[i][defined_indices, :Content]
    end

    return combineData
end

"""
    moving_average_histogram(depth::Vector, Xbar::Vector, windowSize::Int)

Calculate moving average of a given vector and plot histogram with vertical line markers at specified positions

# Arguments:
- `depth`: Input data vector
- `Xbar`: Vector of positions to mark on histogram
- `windowSize`: Moving average window size

# Returns:
- `p`: Plotted figure object
- `moving_avg`: Calculated moving average vector

# Exceptions:
- Throws error when window size exceeds vector length or is less than 1
"""
function moving_average_histogram(depth::Vector, Xbar::Vector, windowSize::Int)
    # Parameter validation
    if windowSize > length(depth)
        throw(ArgumentError("Window size cannot exceed vector length"))
    end
    
    if windowSize < 1
        throw(ArgumentError("Window size must be a positive integer"))
    end
    
    # Calculate moving average
    n = length(depth)
    moving_avg = Vector{Float64}(undef, n - windowSize + 1)
    
    for i in 1:(n - windowSize + 1)
        window = depth[i:(i + windowSize - 1)]
        moving_avg[i] = mean(window)
    end
    
    # Create histogram
    p = histogram(moving_avg, 
                  alpha=0.7,
                  fillcolor=:lightblue)
    
    # Get current plot area y-axis range
    ylims_current = ylims()
    max_height = ylims_current[2]
    
    # Draw red vertical lines at Xbar positions with current window height
    for x_val in Xbar
        # Draw vertical line from y-axis bottom to top
        plot!([x_val, x_val], [0, max_height], 
              color=:red, 
              linewidth=2,
              legend = false)
    end
    
    # Display figure
    display(p)
    
    # Return results
    #return p, moving_avg
end

"""
    paramFineTuning(param, reference_esti, target_esti, refComponent_idx, targetCompnent_idx)

Manually fine-tune parameters to correct inaccurate parameter estimates caused by sample variation

Adjusts means and standard deviations without changing component weights
Adjustment method: translation and scaling based on reference values

# Arguments:
- `param`: Distribution parameters from function fitting_gmm()
- `reference_esti`: Integer indicating position of reference data in param
- `target_esti`: Integer indicating position of data to be corrected in param
- `refComponent_idx`: Vector of reference component positions, must have at least 2 elements
- `targetCompnent_idx`: Vector of target component positions to be corrected

# Returns:
- `param`: Fine-tuned parameters

# Exceptions:
- Throws ArgumentError if input parameter types are incorrect
- Throws ArgumentError if refComponent_idx has fewer than 2 elements
"""
function paramFineTuning(param, reference_esti, target_esti, refComponent_idx, targetCompnent_idx)
    try
        # Input validation
        if !isa(param, Vector)
            throw(ArgumentError("param must be Vector type"))
        end
        
        if !isa(reference_esti, Integer)
            throw(ArgumentError("reference_esti must be Integer type"))
        end
        
        if !isa(target_esti, Integer)
            throw(ArgumentError("target_esti must be Integer type"))
        end
        
        if !isa(refComponent_idx, Vector) || length(refComponent_idx) < 2
            throw(ArgumentError("refComponent_idx must be a vector with at least 2 elements"))
        end
        
                
        # Linear least squares to calculate scaling and translation parameters
        # Using backslash operator to compute y = a * x + b, where a is scaling parameter, b is translation parameter
        x = param[reference_esti].means[refComponent_idx]
        y = param[target_esti].means[refComponent_idx]
        
        num_means = length(refComponent_idx)
        X = hcat(x, ones(num_means))
        coefficients = X \ y  # a = coefficients[1], b = coefficients[2]
        
        # Calculate new values for means and stds, replacing old values
        param[target_esti].means[targetCompnent_idx] = 
            coefficients[1] * param[reference_esti].means[targetCompnent_idx] .+ coefficients[2]
            
        param[target_esti].stds[targetCompnent_idx] = 
            coefficients[1] * param[reference_esti].stds[targetCompnent_idx] .+ coefficients[2]
        
        return param
        
    catch e
        @error "Error occurred in paramFineTuning function: " * string(e)
        rethrow(e)
    end
end

"""
    plotCNVregionMulti(metaDataDF, CNVscaffold, pairedCNV, sample_pair; layout = (1,1), positive_cnv_content = 0.8, chr = "all")

Plot CNV regions for multiple sample pairs

# Arguments:
- `metaDataDF`: Filtered and normalized data frame from function filtering_normalizing_seqDepth()
- `CNVscaffold`: CNV fragment scaffold matrix
- `pairedCNV`: CNV identification for paired samples, 3D array from function copyNumber_comparison()
- `sample_pair`: n×2 matrix representing sample pairs, e.g., [1,2] represents first sample pair
- `layout`: Layout parameter, default is (1, 1)
- `plotsize`: Figure size, default is (1000, 500)
- `positive_cnv_content`: Minimum coverage ratio of genome fragments identified as CNV within CNV region, default is 0.8
- `chr`: Chromosome to plot, default "all" for all chromosomes

# Returns:
- DataFrame containing all plottable CNV regions across all chromosomes
- Figure object that can be displayed directly or further manipulated

# Dependencies:
- Plots: Plotting library
"""
function plotCNVregionMulti(
    metaDataDF,
    CNVscaffold,
    pairedCNV,
    sample_pair;
    layout = (1, 1),
    xlims_pre = [0 0],
    ylims_pre = [0 0],
    plotsize = (1000, 500),
    positive_cnv_content = 0.8,
    chr = "all"
)
    try
        # Input validation
        if !isa(metaDataDF, DataFrame)
            throw(ArgumentError("metaDataDF must be DataFrame type"))
        end
     
        if !isa(CNVscaffold, Matrix) || size(CNVscaffold, 2) != 3
            throw(ArgumentError("CNVscaffold must be an n×3 matrix"))
        end
     
        if !isa(pairedCNV, Array) || length(size(pairedCNV)) != 3
            throw(ArgumentError("pairedCNV must be a 3D array"))
        end
        if !isa(sample_pair, Matrix) || size(sample_pair, 2) != 2
            throw(ArgumentError("sample_pair must be an n×2 matrix"))
        end
     
        num_pairs = size(sample_pair, 1)
        # Create plot
        plt = plot(layout = layout, size = plotsize)
        xlims_current = zeros(num_pairs, 2)
        ylims_current = zeros(num_pairs, 2)
        plotcnv_regionDF = DataFrame()
        for pair_i = 1:num_pairs
            # Data preparation
            cnv_regionsDF_samplepair = cnv_region(
                metaDataDF,
                CNVscaffold,
                pairedCNV[:, sample_pair[pair_i, 1], sample_pair[pair_i, 2]],
                positive_cnv_content
            )
         
            metaDataDF_samplepair = select(
                metaDataDF,
                1, 2, 3,
                3 + sample_pair[pair_i, 1],
                3 + sample_pair[pair_i, 2]
            )
         
            depths_metaData = Matrix(metaDataDF_samplepair[:, 4:end])
            # Calculate mean depth of CNV scaffold fragments
            num_fragments = size(CNVscaffold, 1)
            depths_mean = zeros(num_fragments, 2)
            fragments_chr = String[]
            fragments_chr_startposi = Int64[]
            fragments_chr_endposi = Int64[]
         
            for i = 1:num_fragments
                push!(fragments_chr, String(metaDataDF_samplepair[CNVscaffold[i, 1], 1]))
                push!(fragments_chr_startposi, Int64(metaDataDF_samplepair[CNVscaffold[i, 1], 2]))
                push!(fragments_chr_endposi, Int64(metaDataDF_samplepair[CNVscaffold[i, 2], 3]))
             
                # Add bounds checking
                start_idx = max(1, CNVscaffold[i, 1])
                end_idx = min(size(depths_metaData, 1), CNVscaffold[i, 2])
                if start_idx <= end_idx
                    depths_mean[i, 1] = mean(depths_metaData[start_idx:end_idx, 1])
                    depths_mean[i, 2] = mean(depths_metaData[start_idx:end_idx, 2])
                else
                    depths_mean[i, 1] = 0
                    depths_mean[i, 2] = 0
                end
            end
          
            # Calculate mean_depth_A, mean_depth_B, depth_diff for each CNV region
            if !isempty(cnv_regionsDF_samplepair)
                mean_depth_A = Float64[]
                mean_depth_B = Float64[]
                depth_diff = Float64[]
                for j in 1:nrow(cnv_regionsDF_samplepair)
                    start_bin = cnv_regionsDF_samplepair[j, :start_pos_metaData]
                    end_bin = cnv_regionsDF_samplepair[j, :end_pos_metaData]
                    if start_bin <= end_bin && start_bin >= 1 && end_bin <= size(depths_metaData, 1)
                        mean_A = mean(depths_metaData[start_bin:end_bin, 1])
                        mean_B = mean(depths_metaData[start_bin:end_bin, 2])
                        push!(mean_depth_A, mean_A)
                        push!(mean_depth_B, mean_B)
                        push!(depth_diff, mean_A - mean_B)
                    else
                        push!(mean_depth_A, 0.0)
                        push!(mean_depth_B, 0.0)
                        push!(depth_diff, 0.0)
                    end
                end
                cnv_regionsDF_samplepair[!, :mean_depth_A] = mean_depth_A
                cnv_regionsDF_samplepair[!, :mean_depth_B] = mean_depth_B
                cnv_regionsDF_samplepair[!, :depth_diff] = depth_diff
            end
         
            # Record all plottable regions for all chromosomes
            sample1 = repeat([sample_pair[pair_i,1]], size(cnv_regionsDF_samplepair, 1))
            sample2 = repeat([sample_pair[pair_i,2]], size(cnv_regionsDF_samplepair, 1))
            temp_DF = DataFrame(pair_sampleA = sample1, pair_sampleB = sample2)
            temp_DF = hcat(temp_DF, cnv_regionsDF_samplepair)
            plotcnv_regionDF = vcat(plotcnv_regionDF, temp_DF)
         
            # Select fragments to plot
            if chr == "all"
                # Remove fragments that don't meet requirements
                selected_idx = findall(fragments_chr_endposi .- fragments_chr_startposi .> 0)
                # Plot all chromosomes
                temp_x = 1:length(selected_idx)
                temp_y1 = depths_mean[selected_idx, 1]
                plot!(
                    plt,
                    subplot = pair_i,
                    temp_x,
                    temp_y1,
                    linewidth = 0.5,
                    color = "red",
                    legend = false,
                    label = "Sample $(sample_pair[pair_i, 1])"
                )
             
                temp_y2 = depths_mean[selected_idx, 2]
                plot!(
                    plt,
                    subplot = pair_i,
                    temp_x,
                    temp_y2,
                    linewidth = 0.5,
                    color = "blue",
                    legend = false,
                    label = "Sample $(sample_pair[pair_i, 2])"
                )
             
                # If there is CNV region data, plot CNV regions
                if size(cnv_regionsDF_samplepair, 1) > 0
                    # Adjust according to actual data structure
                    # Example: Draw black lines to represent CNV regions
                    # plot!(plt, subplot=pair_i, temp_xCNV', temp_yCNV', linewidth = 1, color = "black", label = "CNV Regions")
                end
            else
                # Plot specified chromosome
                selected_idx = findall((fragments_chr .== chr) .&& (fragments_chr_endposi .- fragments_chr_startposi .> 0))
                if !isempty(selected_idx)
                    temp_x = hcat(fragments_chr_startposi[selected_idx], fragments_chr_endposi[selected_idx])
                    temp_y1 = hcat(depths_mean[selected_idx, 1], depths_mean[selected_idx, 1])
                    plot!(
                        plt,
                        subplot = pair_i,
                        temp_x',
                        temp_y1',
                        linewidth = 1,
                        color = "red",
                        legend = false,
                        title = "Sample $(sample_pair[pair_i, :])"
                    )
                 
                    temp_y2 = hcat(depths_mean[selected_idx, 2], depths_mean[selected_idx, 2])
                    plot!(
                        plt,
                        subplot = pair_i,
                        temp_x',
                        temp_y2',
                        linewidth = 1,
                        color = "blue",
                        legend = false,
                    )
                 
                    # Plot CNV regions
                    if isa(cnv_regionsDF_samplepair, DataFrame) && size(cnv_regionsDF_samplepair, 1) > 0 && size(cnv_regionsDF_samplepair, 2) >= 3
                        chr_cnv_rows = cnv_regionsDF_samplepair[:, :chr] .== chr
                        if any(chr_cnv_rows)
                            selected_cnv_idx = findall(chr_cnv_rows)
                            # Filter out rows where |depth_diff| < 0
                            valid_mask = abs.(cnv_regionsDF_samplepair[selected_cnv_idx, :depth_diff]) .>= 0
                            selected_cnv_idx = selected_cnv_idx[valid_mask]
                            if !isempty(selected_cnv_idx)
                                temp_xCNV = hcat(
                                    cnv_regionsDF_samplepair[selected_cnv_idx, :start_pos],
                                    cnv_regionsDF_samplepair[selected_cnv_idx, :end_pos]
                                )
                                temp_yCNV = ones(length(selected_cnv_idx), 2) * 1 # Fixed Y value for visualization
                                plot!(
                                    plt,
                                    subplot = pair_i,
                                    temp_xCNV',
                                    temp_yCNV',
                                    linewidth = 2,
                                    color = "black",
                                    legend = false
                                )
                            end
                        end
                    end
                else
                    @warn "No data found for chromosome $chr"
                end
            end
        end
     
        # Adjust display range
        for pair_i = 1:num_pairs # Get current parameters
            xlims_current[pair_i, :] = [xlims(plt[pair_i])[1] xlims(plt[pair_i])[2]]
            ylims_current[pair_i, :] = [ylims(plt[pair_i])[1] ylims(plt[pair_i])[2]]
        end
        min_xlims = minimum(xlims_current[:, 1]) # Get extreme values
        min_ylims = minimum(ylims_current[:, 1])
        max_ylims = maximum(ylims_current[:, 2])
        max_xlims = maximum(xlims_current[:, 2])
     
        if xlims_pre == [0 0]
            for pair_i = 1:num_pairs # Adjust view
                xlims!(plt[pair_i], min_xlims, max_xlims)
            end
        else
            for pair_i = 1:num_pairs # Adjust view
                xlims!(plt[pair_i], xlims_pre[1], xlims_pre[2])
            end
        end
        if ylims_pre == [0 0]
            for pair_i = 1:num_pairs # Adjust view
                ylims!(plt[pair_i], min_ylims, max_ylims)
            end
        else
            for pair_i = 1:num_pairs # Adjust view
                ylims!(plt[pair_i], ylims_pre[1], ylims_pre[2])
            end
        end
         
        # Display figure
        display(plt)
       
        # Create filtered subplotcnv_regionDF
        if hasproperty(plotcnv_regionDF, :depth_diff) && !isempty(plotcnv_regionDF)
            valid_filter = abs.(coalesce.(plotcnv_regionDF.depth_diff, 0.0)) .>= 20
            subplotcnv_regionDF = plotcnv_regionDF[valid_filter, :]
        else
            subplotcnv_regionDF = DataFrame()
        end
       
        return plotcnv_regionDF, subplotcnv_regionDF, plt
     
    catch e
        @error "Error occurred in plotCNVregionMulti function: " * string(e)
        rethrow(e)
    end
end








"""
    plotCNVregion(metaDataDF, CNVscaffold, pairedCNV, sample_pair; positive_cnv_content = 0.8, chr = "all")

Plot CNV region diagram

# Arguments:
- `metaDataDF`: Filtered and normalized data frame from function filtering_normalizing_seqDepth()
- `CNVscarfold`: CNV fragment scaffold matrix
- `pairedCNV`: CNV identification for paired samples, 3D array from function copyNumber_comparison()
- `sample_pair`: Vector with two numeric elements representing sample pair, e.g., [1,2] represents first and second sample pair
- `positive_cnv_content`: Minimum coverage ratio of genome fragments identified as CNV within CNV region, default is 0.8
- `chr`: Chromosome to plot, default "all" for all chromosomes

# Returns:
- No return value, displays figure directly

# Dependencies:
- Plots: Plotting library
"""
function plotCNVregion(metaDataDF, CNVscarfold, pairedCNV, sample_pair; positive_cnv_content = 0.8, chr = "all")
    
    try
        # Input validation
        if !isa(metaDataDF, DataFrame)
            throw(ArgumentError("metaDataDF must be DataFrame type"))
        end
        
        if !isa(CNVscarfold, Matrix) || size(CNVscarfold, 2) != 3
            throw(ArgumentError("CNVscarfold must be an n×3 matrix"))
        end               
        
        if !isa(pairedCNV, Array) || length(size(pairedCNV)) != 3
            throw(ArgumentError("pairedCNV must be a 3D array"))
        end

        if !isa(sample_pair, Vector) || length(sample_pair) != 2
            throw(ArgumentError("sample_pair must be a vector of length 2"))
        end  
        

        # Data preparation        
        cnv_regionsDF_samplepair = cnv_region(metaDataDF, CNVscarfold, pairedCNV[:, sample_pair[1], sample_pair[2]], positive_cnv_content)
        metaDataDF_samplepair = select(metaDataDF, 1, 2, 3, 3 + sample_pair[1], 3 + sample_pair[2])        
        depths_metaData = Matrix(metaDataDF_samplepair[:, 4:end])

        # Calculate mean depth of CNV scaffold fragments
        num_fragments = size(CNVscarfold, 1)
        depths_mean = zeros(num_fragments, 2)
        fragments_chr = String[]
        fragments_chr_startposi = Int64[]
        fragments_chr_endposi = Int64[]
        
        for i = 1:num_fragments
            push!(fragments_chr, String(metaDataDF_samplepair[CNVscarfold[i, 1], 1]))
            push!(fragments_chr_startposi, Int64(metaDataDF_samplepair[CNVscarfold[i, 1], 2]))
            push!(fragments_chr_endposi, Int64(metaDataDF_samplepair[CNVscarfold[i, 2], 3]))
            # Add bounds checking
            start_idx = max(1, CNVscarfold[i, 1])
            end_idx = min(size(depths_metaData, 1), CNVscarfold[i, 2])
            if start_idx <= end_idx
                depths_mean[i, 1] = mean(depths_metaData[start_idx:end_idx, 1])
                depths_mean[i, 2] = mean(depths_metaData[start_idx:end_idx, 2])
            else
                depths_mean[i, 1] = 0
                depths_mean[i, 2] = 0
            end
        end

        # Create plot
        plt = plot()
        
        # Select fragments to plot
        if chr == "all"
            # Remove fragments that don't meet requirements
            selected_idx = findall(fragments_chr_endposi .- fragments_chr_startposi .> 0)

            # Plot all chromosomes
            temp_x = 1:length(selected_idx)
            temp_y1 = depths_mean[selected_idx, 1]
            plot!(plt, temp_x, temp_y1, linewidth = 0.5, color = "red", legend = false)
            temp_y2 = depths_mean[selected_idx, 2]
            plot!(plt, temp_x, temp_y2, linewidth = 0.5, color = "blue", legend = false)
            
            # If there is CNV region data, plot CNV regions
            if size(cnv_regionsDF_samplepair, 1) > 0
                # Adjust according to actual data structure
                # Example: Draw red lines to represent CNV regions
                # plot!(plt, temp_xCNV', temp_yCNV', linewidth = 1, color = "red", label = "CNV Regions")
            end
        else
            # Plot specified chromosome
            selected_idx = findall((fragments_chr .== chr) .&& (fragments_chr_endposi .- fragments_chr_startposi .> 0))
            if !isempty(selected_idx)
                temp_x = hcat(fragments_chr_startposi[selected_idx], fragments_chr_endposi[selected_idx])
                temp_y1 = hcat(depths_mean[selected_idx, 1], depths_mean[selected_idx, 1])
                plot!(plt, temp_x', temp_y1', linewidth = 1, color = "red", legend = false)
                temp_y2 = hcat(depths_mean[selected_idx, 2], depths_mean[selected_idx, 2])
                plot!(plt, temp_x', temp_y2', linewidth = 1, color = "blue", legend = false)
                
                # Plot CNV regions
                chr_cnv_rows = cnv_regionsDF_samplepair[:, 1] .== chr
                if any(chr_cnv_rows) && size(cnv_regionsDF_samplepair, 2) >= 3
                    selected_cnv_idx = findall(chr_cnv_rows)
                    temp_xCNV = hcat(
                        cnv_regionsDF_samplepair[selected_cnv_idx, 2],
                        cnv_regionsDF_samplepair[selected_cnv_idx, 3]
                    )
                    temp_yCNV = ones(length(selected_cnv_idx), 2) * 1  # Fixed Y value for visualization
                    plot!(plt, temp_xCNV', temp_yCNV', linewidth = 2, color = "black", legend = false )
                end
            else
                @warn "No data found for chromosome $chr"
            end
        end
        
        # Display figure
        display(plt)
        return plt
        
    catch e
        @error "Error occurred in plotCNVregion function: " * string(e)
        rethrow(e)
    end
end


"""
    cnv_region(metaDataDF, CNVscarffold, fragment_status, threshold_proportion)

Synthesize CNV regions based on fragment status

# Arguments:
- `metaDataDF`: Filtered and normalized data frame from function filtering_normalizing_seqDepth()
- `CNVscarffold`: Fragment information matrix from function fragment_status_matrix()
- `fragment_status`: Fragment CNV status, elements correspond to rows of CNVscarffold, can be from copyNumber_comparison
- `threshold_proportion`: Fragment CNV content threshold, e.g., 0.8

# Returns:
- `coordinate_cnv_regions`: CNV region coordinate information

# Dependencies:
- synth_cnv: CNV synthesis function

# Exceptions:
- Throws ArgumentError if input parameter types are incorrect
- Throws DimensionMismatch if parameter dimensions don't match
"""
function cnv_region(metaDataDF, CNVscarffold, fragment_status, threshold_proportion)
    try
        # Input validation
        if !isa(metaDataDF, DataFrame)
            throw(ArgumentError("metaDataDF must be DataFrame type"))
        end
        
        if !isa(CNVscarffold, Matrix) || size(CNVscarffold, 2) != 3
            throw(ArgumentError("CNVscarffold must be an n×3 matrix"))
        end
        
        if !isa(fragment_status, Vector)
            throw(ArgumentError("fragment_status must be Vector type"))
        end
        
        if !isa(threshold_proportion, Number) || threshold_proportion < 0 || threshold_proportion > 1
            throw(ArgumentError("threshold_proportion must be a number in [0,1] range"))
        end
        
        # Check DataFrame has enough columns (at least 3: chr, binLowEdge, binUpEdge)
        if size(metaDataDF, 2) < 3
            throw(DimensionMismatch("metaDataDF requires at least 3 columns"))
        end
        
        # Check dimension compatibility
        if length(fragment_status) != size(CNVscarffold, 1)
            throw(DimensionMismatch("fragment_status length must match CNVscarffold row count"))
        end
        
        # Build fragment coordinate matrix
        # Column 1: chromosome ID, Column 2: start position, Column 3: end position, Column 4: fragment length
        fragment_coords = hcat(
            metaDataDF[CNVscarffold[:, 1], 1],  # Chromosome ID
            CNVscarffold[:, 1],                 # Start position
            CNVscarffold[:, 2],                 # End position
            CNVscarffold[:, 2] - CNVscarffold[:, 1] .+ 1  # Fragment length
        )
        
        # Use synth_cnv function to synthesize CNV regions
        temp_regions = synth_cnv(fragment_coords, fragment_status, threshold_proportion)
        
        # If synth_cnv returns error code, return directly
        if temp_regions == -1
            @error "Error occurred during CNV synthesis"
            return []
        end
        

        # Extract CNV region coordinate information
        chr = String[]
        start_pos = Int64[]
        end_pos = Int64[]
        num_regions = size(temp_regions, 1)
        start_pos_metaData = Int64[]
        end_pos_metaData = Int64[]         
        if !isempty(temp_regions)
            for i in 1:num_regions
                temp_chr = String(metaDataDF[CNVscarffold[temp_regions[i, 1], 1], 1])
                temp_start_pos = Int64(metaDataDF[CNVscarffold[temp_regions[i, 1], 1], 2])
                temp_end_pos = Int64(metaDataDF[CNVscarffold[temp_regions[i, 2], 2], 3])
                if temp_end_pos >= temp_start_pos # Ensure indices don't exceed bounds
                    push!(chr, temp_chr)
                    push!(start_pos, temp_start_pos)
                    push!(end_pos, temp_end_pos)
                    push!(start_pos_metaData, CNVscarffold[temp_regions[i, 1], 1])
                    push!(end_pos_metaData, CNVscarffold[temp_regions[i, 2], 2])
                end
            end
        end

        coordinate_cnv_regions = DataFrame()
        coordinate_cnv_regions[!,:chr] = chr
        coordinate_cnv_regions[!,:start_pos] = start_pos
        coordinate_cnv_regions[!,:end_pos] = end_pos
        coordinate_cnv_regions[!,:length] = end_pos .- start_pos
        coordinate_cnv_regions[!,:start_pos_metaData] = start_pos_metaData
        coordinate_cnv_regions[!,:end_pos_metaData] = end_pos_metaData

        return coordinate_cnv_regions
        
    catch e
        @error "Error occurred in cnv_region function: " * string(e)
        rethrow(e)
    end
end


"""
    copyNumber_comparison(metaDataDF, CNVscaffold, copynumber_info, gmmParam, windowSize; pvalue_theSame = 0.05, fdr_threshold = 0.2)

Identify copy number differences between samples.

# Arguments:
- `metaDataDF`: Filtered and normalized data frame, output from `filtering_normalizing_seqDepth()`
- `CNVscaffold`: CNV candidate interval matrix, generated by `genome_segment()`
- `copynumber_info`: Copy number estimation matrix, generated by `copynumber_determination()`
- `gmmParam`: GMM parameter vector, generated by `fitting_gmm()`
- `windowSize`: Window size used when estimating GMM parameters
- `pvalue_theSame`: P-value threshold for copy number consistency test (default: 0.05)
- `fdr_threshold`: FDR threshold, lenient to reduce false negatives (default: 0.2)

# Returns:
- `cnv_state_info`: 3D array indicating CNV state differences between samples
    - Dimension 1: CNV fragments
    - Dimension 2: Sample j
    - Dimension 3: Sample k

# Dependencies:
- `prob_copy_differenceM`: Custom function to compute probability of copy number difference
- `upgrade_state_vector`: Custom function to update state vector

# Exception handling:
- Throws `ArgumentError` if input types are incorrect
- Throws `DimensionMismatch` if dimensions do not match
"""
function copyNumber_comparison(metaDataDF, CNVscaffold, copynumber_info, gmmParam, windowSize; pvalue_theSame = 0.05, fdr_threshold = 0.1)
    try
        # Input validation
        if !isa(metaDataDF, DataFrame)
            throw(ArgumentError("metaDataDF must be a DataFrame"))
        end
        
        if !isa(CNVscaffold, Matrix) || size(CNVscaffold, 2) != 3
            throw(ArgumentError("CNVscaffold must be a matrix with n rows and 3 columns"))
        end
        
        if !isa(copynumber_info, Array) || ndims(copynumber_info) != 3
            throw(ArgumentError("copynumber_info must be a 3-dimensional array"))
        end
        
        if !isa(gmmParam, Vector)
            throw(ArgumentError("gmmParam must be a Vector"))
        end
        
        if !isa(pvalue_theSame, Number) || pvalue_theSame <= 0 || pvalue_theSame >= 1
            throw(ArgumentError("pvalue_theSame must be a number in (0,1)"))
        end
        
        if !isa(fdr_threshold, Number) || fdr_threshold <= 0 || fdr_threshold >= 1
            throw(ArgumentError("fdr_threshold must be a number in (0,1)"))
        end
        
        # Check if DataFrame has enough columns (at least 4: chr, binLowEdge, binUpEdge, depth data)
        if size(metaDataDF, 2) < 4
            throw(DimensionMismatch("metaDataDF must have at least 4 columns"))
        end
        
        # Check dimension compatibility
        num_CNVfragment = size(CNVscaffold, 1)
        num_sample = size(copynumber_info, 1)
        
        if size(copynumber_info, 2) != num_CNVfragment
            throw(DimensionMismatch("Second dimension of copynumber_info must match number of rows in CNVscaffold"))
        end
        
        if length(gmmParam) != num_sample
            throw(DimensionMismatch("Length of gmmParam must match the number of samples"))
        end
        
        # **************** Identify copy number differences between samples *****************
        subData = Matrix(metaDataDF[:, 4:end])
        num_fragment = size(subData, 1)
        num_samples = size(subData, 2)
        
        # Check if indices in CNVscaffold are within valid range
        if !isempty(CNVscaffold)
            max_index = maximum(CNVscaffold[:, [1, 2]])
            if max_index > num_fragment
                throw(BoundsError("Indices in CNVscaffold exceed the number of rows in subData"))
            end
            
            min_index = minimum(CNVscaffold[:, [1, 2]])
            if min_index < 1
                throw(BoundsError("Indices in CNVscaffold cannot be less than 1"))
            end
        end
        
        cnv_state_info = zeros(num_CNVfragment, num_sample, num_sample)  # 3D array for state
        cnv_prob_info  = zeros(num_CNVfragment, num_sample, num_sample)  # 3D array for probabilities
        
        for j in 1:num_sample
            for k in 1:num_sample
                for i in 1:num_CNVfragment
                    # Both fragments are confidently called and have different copy numbers
                    if (copynumber_info[j, i, 1] > 0) && (copynumber_info[k, i, 1] > 0) &&
                       (copynumber_info[j, i, 1] != copynumber_info[k, i, 1])
                        try
                            temp_star = CNVscaffold[i, 1]
                            temp_end  = CNVscaffold[i, 2]
                            # Ensure indices are within bounds
                            frag_start = max(1, temp_star)
                            frag_end   = min(num_fragment, temp_end)
                            
                            # Calculate probability of copy number difference
                            cnv_prob_info[i, j, k] = prob_copy_differenceM(
                                [subData[frag_start:frag_end, j] subData[frag_start:frag_end, k]],
                                [copynumber_info[j, i, 1] copynumber_info[k, i, 1]],
                                gmmParam[[j k]],
                                windowSize
                            )
                            
                            # Mark as different if probability of same copy number is low
                            if (cnv_prob_info[i, j, k] < pvalue_theSame) &&
                               (copynumber_info[j, i, 1] != copynumber_info[k, i, 1])
                                cnv_state_info[i, j, k] = 1
                            end
                        catch e
                            @warn "Error calculating copy number difference for fragment $i between samples $j and $k: " * string(e)
                        end
                    else
                        cnv_prob_info[i, j, k] = 1.0
                    end
                end
                # Optional debug prints (commented out)
                # println("Copy number difference comparison completed for samples $j and $k: ", sum(cnv_state_info[:, j, k]))
                # println("Number of zero probabilities for samples $j and $k: ", sum(cnv_prob_info[:, j, k] .== 0))
            end
        end
        
        # ************* Expand CNV difference regions (lenient FDR to reduce false negatives) *************
        # The following block is currently disabled / under improvement - comment preserved for reference
        
        for j in 1:num_sample
            for k in 1:num_sample
                try
                    while true
                        temp_state_info = upgrade_state_vector(cnv_prob_info[:, j, k], cnv_state_info[:, j, k], fdr_threshold)
                        
                        if sum(temp_state_info .!= cnv_state_info[:, j, k]) > 0
                            cnv_state_info[:, j, k] = temp_state_info
                        else
                            break
                        end
                    end
                catch e
                    @warn "Error updating state vector for samples $j and $k: " * string(e)
                end
            end
        end
        
        return cnv_state_info  # , cnv_prob_info  (commented out - return only state by default)
        
    catch e
        @error "Error occurred in copyNumber_comparison function: " * string(e)
        rethrow(e)
    end
end



function copyNumber_comparison_filter(metaDataDF, CNVscaffold, copynumber_info, gmmParam, windowSize; pvalue_theSame = 0.05, fdr_threshold = 0.1 , highRDthreshold = 100)
    try
        # Input validation
        if !isa(metaDataDF, DataFrame)
            throw(ArgumentError("metaDataDF must be DataFrame type"))
        end
        
        if !isa(CNVscaffold, Matrix) || size(CNVscaffold, 2) != 3
            throw(ArgumentError("CNVscaffold must be an n×3 matrix"))
        end
        
        if !isa(copynumber_info, Array) || ndims(copynumber_info) != 3
            throw(ArgumentError("copynumber_info must be a 3D array"))
        end
        
        if !isa(gmmParam, Vector)
            throw(ArgumentError("gmmParam must be Vector type"))
        end
        
        if !isa(pvalue_theSame, Number) || pvalue_theSame <= 0 || pvalue_theSame >= 1
            throw(ArgumentError("pvalue_theSame must be a number in (0,1) range"))
        end
        
        if !isa(fdr_threshold, Number) || fdr_threshold <= 0 || fdr_threshold >= 1
            throw(ArgumentError("fdr_threshold must be a number in (0,1) range"))
        end
        
        # Check DataFrame has enough columns (at least 4: chr, binLowEdge, binUpEdge, depth data)
        if size(metaDataDF, 2) < 4
            throw(DimensionMismatch("metaDataDF requires at least 4 columns"))
        end
        
        # Check dimension compatibility
        num_CNVfragment = size(CNVscaffold, 1)
        num_sample = size(copynumber_info, 1)
        
        if size(copynumber_info, 2) != num_CNVfragment
            throw(DimensionMismatch("copynumber_info dimension 2 must match CNVscaffold row count"))
        end
        
        if length(gmmParam) != num_sample
            throw(DimensionMismatch("gmmParam length must match sample count"))
        end
        
        # ****************Identify copy number differences between samples*****************
        subData = Matrix(metaDataDF[:, 4:end])
        num_fragment = size(subData, 1)
        num_samples = size(subData, 2)
        
        # Check CNVscaffold indices are within valid range
        if !isempty(CNVscaffold)
            max_index = maximum(CNVscaffold[:, [1, 2]])
            if max_index > num_fragment
                throw(BoundsError("CNVscaffold indices exceed subData row count"))
            end
            
            min_index = minimum(CNVscaffold[:, [1, 2]])
            if min_index < 1
                throw(BoundsError("CNVscaffold indices cannot be less than 1"))
            end
        end
    
        
        cnv_state_info = zeros(num_CNVfragment, num_sample, num_sample)  # 3D table
        cnv_prob_info = zeros(num_CNVfragment, num_sample, num_sample)   # 3D table
        
        
        for j in 1:num_sample
            for k in 1:num_sample
                for i in 1:num_CNVfragment                    
                    # Both fragments correctly identified with same copy number (retain probability judgment only, don't use deltaRDmean threshold directly)
                    if (copynumber_info[j, i, 1] > 0) & (copynumber_info[k, i, 1] > 0) & 
                       (copynumber_info[j, i, 1] == copynumber_info[k, i, 1])
                        # This branch no longer directly sets cnv_state_info=1 (moved to extension phase)
                        continue
                    elseif (copynumber_info[j, i, 1] > 0) & (copynumber_info[k, i, 1] > 0) & 
                       (copynumber_info[j, i, 1] != copynumber_info[k, i, 1])
                        try
                            temp_star = CNVscaffold[i, 1]
                            temp_end = CNVscaffold[i, 2]
                            # Ensure indices are within valid range
                            frag_start = max(1, temp_star)
                            frag_end = min(num_fragment, temp_end)

                            # Calculate copy number difference probability (core judgment condition)
                            cnv_prob_info[i, j, k] = prob_copy_differenceM(
                                [subData[frag_start:frag_end, j] subData[frag_start:frag_end, k]], 
                                [copynumber_info[j, i, 1] copynumber_info[k, i, 1]], 
                                gmmParam[[j k]],
                                windowSize
                            )
                                                       
                            if (cnv_prob_info[i, j, k] < pvalue_theSame)
                                cnv_state_info[i, j, k] = 1
                            end 
                        catch e
                            @warn "Error calculating copy number difference for fragment $i samples $j and $k: " * string(e)
                        end
                    elseif (copynumber_info[j, i, 1] < 0) | (copynumber_info[k, i, 1] < 0)
                        # This branch no longer directly sets cnv_state_info=1 (moved to extension phase)
                        continue
                    else
                         cnv_prob_info[i, j, k] = 1
                    end                    
                end
            end
        end
        
        # *************Extend copy number difference regions (relax false positives, reduce false negatives)*****************
        # First obtain initial cnv_state_info based on core probability judgment
        # Then extend: for fragments originally 0, if deltaRDmean > highRDthreshold and at least one adjacent fragment is 1, extend to 1
        # Iterate multiple times here (3 times is enough to propagate most contiguous regions), avoid incomplete propagation from single traversal
        for iter in 1:3
            for j in 1:num_sample
                for k in 1:num_sample
                    for i in 1:num_CNVfragment
                        if cnv_state_info[i, j, k] == 0  # Only extend those originally 0
                            temp_star = CNVscaffold[i, 1]
                            temp_end = CNVscaffold[i, 2]
                            frag_start = max(1, temp_star)
                            frag_end = min(num_fragment, temp_end)
                            
                            deltaRDmean = abs(mean(subData[frag_start:frag_end, j]) - mean(subData[frag_start:frag_end, k]))
                            
                            if deltaRDmean > highRDthreshold
                                has_neighbor_diff = false
                                if i > 1 && cnv_state_info[i-1, j, k] == 1
                                    has_neighbor_diff = true
                                end
                                if i < num_CNVfragment && cnv_state_info[i+1, j, k] == 1
                                    has_neighbor_diff = true
                                end
                                
                                if has_neighbor_diff
                                    cnv_state_info[i, j, k] = 1
                                end
                            end
                        end
                    end
                end
            end
        end
                
        return cnv_state_info #, cnv_prob_info
        
    catch e
        @error "Error occurred in copyNumber_comparison function: " * string(e)
        rethrow(e)
    end
end


# First perform copy number difference identification, then run the following function to reduce RD false negatives
# Usage: pairedCNVT293FN = copyNumber_comparison_filter_fn(metaDataDF, T293CNVscaffold, copynumber_info, pairedCNVT293; highRDthreshold=20);
# pairedCNVT293FN is a supplemented version of initially identified paired sample copy number difference data pairedCNVT293
# Added region R: adjacent to positive region, with RD difference greater than threshold (includes all regardless of copy number identification being same, different, or -1)
# Iterate multiple rounds until no further extension possible
function copyNumber_comparison_filter_fn(metaDataDF, CNVscaffold, copynumber_info, cnv_state_info;   highRDthreshold = 100)
    try
        # Input validation
        if !isa(metaDataDF, DataFrame)
            throw(ArgumentError("metaDataDF must be DataFrame type"))
        end
        
        if !isa(CNVscaffold, Matrix) || size(CNVscaffold, 2) != 3
            throw(ArgumentError("CNVscaffold must be an n×3 matrix"))
        end
        
        if !isa(copynumber_info, Array) || ndims(copynumber_info) != 3
            throw(ArgumentError("copynumber_info must be a 3D array"))
        end
       
        if !isa(cnv_state_info, Array) || ndims(copynumber_info) != 3
            throw(ArgumentError("copynumber_info must be a 3D array"))
        end
        
        
        # Check DataFrame has enough columns (at least 4: chr, binLowEdge, binUpEdge, depth data)
        if size(metaDataDF, 2) < 4
            throw(DimensionMismatch("metaDataDF requires at least 4 columns"))
        end
        
        # Check dimension compatibility
        num_CNVfragment = size(CNVscaffold, 1)
        num_sample = size(copynumber_info, 1)
        
        if size(copynumber_info, 2) != num_CNVfragment
            throw(DimensionMismatch("copynumber_info dimension 2 must match CNVscaffold row count"))
        end
        
    
        
        # ****************Identify copy number differences between samples*****************
        subData = Matrix(metaDataDF[:, 4:end])
        num_fragment = size(subData, 1)
        num_samples = size(subData, 2)
        
        # Check CNVscaffold indices are within valid range
        if !isempty(CNVscaffold)
            max_index = maximum(CNVscaffold[:, [1, 2]])
            if max_index > num_fragment
                throw(BoundsError("CNVscaffold indices exceed subData row count"))
            end
            
            min_index = minimum(CNVscaffold[:, [1, 2]])
            if min_index < 1
                throw(BoundsError("CNVscaffold indices cannot be less than 1"))
            end
        end
    
        
        # Initialize: first round input is original cnv_state_info
        cnv_state_info_prev = copy(cnv_state_info)
        # Flag for whether to continue loop (initially true, execute at least one round)
        has_change = true
        # Record iteration count (optional, for monitoring)
        iter_num = 0

        # Iterative loop: continue as long as there are updates
        while has_change
            iter_num += 1
            @info "Starting iteration round $iter_num"
            # Initialize new state for this round (based on previous round results)
            cnv_state_info_new = copy(cnv_state_info_prev)
            # Flag for whether this round has updates
            current_change = false
            # Core logic loop (same as original code, just replace cnv_state_info with cnv_state_info_prev)
            for j in 1:num_samples
                for k in 1:num_samples 
                    for i in 2:(num_CNVfragment-1)
                        # Both fragments correctly identified with same copy number
                        if (copynumber_info[j, i, 1] > 0) & (copynumber_info[k, i, 1] > 0)
                            try
                                temp_star = CNVscaffold[i, 1]
                                temp_end = CNVscaffold[i, 2]
                                # Ensure indices are within valid range
                                frag_start = max(1, temp_star)
                                frag_end = min(num_fragment, temp_end)
            
                                deltaRDmean = abs(mean(subData[frag_start:frag_end, j]) - mean(subData[frag_start:frag_end, k]))
                                if (deltaRDmean .> highRDthreshold) & ((cnv_state_info_prev[i-1, j, k] == 1) | (cnv_state_info_prev[i+1, j, k] == 1))
                                # Check if different from original value, if so mark as updated
                                   if (cnv_state_info_new[i, j, k] != 1)
                                       cnv_state_info_new[i, j, k] = 1
                                       current_change = true
                                   end
                                end 
                            catch e
                                @warn "Error calculating copy number difference for fragment $i samples $j and $k: " * string(e)
                            end
                
                        elseif (copynumber_info[j, i, 1] < 0) | (copynumber_info[k, i, 1] < 0)
                            try
                                temp_star = CNVscaffold[i, 1]
                                temp_end = CNVscaffold[i, 2]
                                # Ensure indices are within valid range
                                frag_start = max(1, temp_star)
                                frag_end = min(num_fragment, temp_end)
                                deltaRDmean = abs(mean(subData[frag_start:frag_end, j]) - mean(subData[frag_start:frag_end, k]))
                                if (deltaRDmean .> highRDthreshold) & 
                                    ((cnv_state_info_prev[i-1, j, k] == 1) | (cnv_state_info_prev[i+1, j, k] == 1))
                                # Check if different from original value, if so mark as updated
                                    if (cnv_state_info_new[i, j, k] != 1)
                                        cnv_state_info_new[i, j, k] = 1
                                        current_change = true
                                    end
                                end
                            catch e
                                @warn "Error calculating copy number difference for fragment $i samples $j and $k: " * string(e)
                            end
                        end                    
                    end
        
                end
            end

    # Update: use this round's result as input for next round
            cnv_state_info_prev = copy(cnv_state_info_new)
    # Update loop flag: continue if this round had updates, otherwise terminate
            has_change = current_change

            @info "Iteration round $iter_num completed, has updates: $has_change"
        end

        # Final result stored in cnv_state_info_prev (or cnv_state_info_new)
        @info "Iteration finished, total iteration rounds: $iter_num"
        final_cnv_state_info = cnv_state_info_prev
        
        # *************Extend copy number difference regions, relax false positives, reduce false negatives*****************
        # The following processing may cause unnecessary complications, temporarily abandoned, to be improved
        
                
        return final_cnv_state_info #, cnv_prob_info
        
    catch e
        @error "Error occurred in copyNumber_comparison function: " * string(e)
        rethrow(e)
    end
end


"""
    copynumber_determination(metaDataDF, CNVscaffold, gmmParam, flanking_length; thresholdBFactor = 2.0)

Identify copy number of genomic fragments

# Arguments:
- `metaDataDF`: Filtered and normalized data frame from function filtering_normalizing_seqDepth()
- `CNVscaffold`: CNV examination interval scaffold matrix, generated by genome_segment function
- `gmmParam`: Model parameters from function fitting_gmm()
- `flanking_length`: Left and right boundary length of fragments, must be consistent with genome_segment() call
- `thresholdBFactor`: Bayes factor threshold, typical default value is 2.0

# Returns:
- `copyNumberInfo`: Copy number information, 3D table, 1D samples, 2D fragments, 3D copy number details
  Copy number details: first column is copy number, subsequent columns are probabilities for copy numbers k-1

# Dependencies:
- identify_copy_number: Custom function to identify copy number

# Exceptions:
- Throws ArgumentError if input parameter types are incorrect
- Throws DimensionMismatch if parameter dimensions don't match
- Throws BoundsError if indices exceed bounds
"""
function copynumber_determination(metaDataDF, CNVscaffold, gmmParam, flanking_length; thresholdBFactor = 2.0)
    try
        # Input validation
        if !isa(metaDataDF, DataFrame)
            throw(ArgumentError("metaDataDF must be DataFrame type"))
        end
        
        if !isa(CNVscaffold, Matrix) || size(CNVscaffold, 2) != 3
            throw(ArgumentError("CNVscaffold must be an n×3 matrix"))
        end
        
        if !isa(gmmParam, Vector)
            throw(ArgumentError("gmmParam must be Vector type"))
        end
        
        if !isa(flanking_length, Integer) || flanking_length <= 0
            throw(ArgumentError("flanking_length must be a positive integer"))
        end
        
        if !isa(thresholdBFactor, Number) || thresholdBFactor <= 0
            throw(ArgumentError("thresholdBFactor must be a positive number"))
        end
        
        # Check DataFrame has enough columns (at least 4: chr, binLowEdge, binUpEdge, depth data)
        if size(metaDataDF, 2) < 4
            throw(DimensionMismatch("metaDataDF requires at least 4 columns"))
        end
        
        # ***************Determine fragment copy number*******************
        subData = Matrix(metaDataDF[:, 4:end])
        num_row = size(subData, 1)
        num_col = size(subData, 2)
        
        # Check parameter compatibility
        if flanking_length * 2 >= num_row
            throw(ArgumentError("flanking_length too large, incompatible with data length"))
        end
        
        if length(gmmParam) != num_col
            throw(DimensionMismatch("gmmParam length must match sample count"))
        end
        
        # Check CNVscaffold indices are within valid range
        if !isempty(CNVscaffold)
            max_index = maximum(CNVscaffold[:, [1, 2]])
            if max_index > num_row
                throw(BoundsError("CNVscaffold indices exceed subData row count"))
            end
            
            min_index = minimum(CNVscaffold[:, [1, 2]])
            if min_index < 1
                throw(BoundsError("CNVscaffold indices cannot be less than 1"))
            end
        end
        
        AllMeanSub = zeros(size(subData))
        for i in (flanking_length + 1):(num_row - flanking_length)
            for j in 1:num_col
                start_idx = max(1, i - flanking_length)
                end_idx = min(num_row, i + flanking_length)
                AllMeanSub[i, j] = mean(subData[start_idx:end_idx, j])
            end
        end
        
        depthFrag = zeros(size(CNVscaffold, 1), num_col)
        for i in 1:size(CNVscaffold, 1)
            for j in 1:num_col
                temp_star = CNVscaffold[i, 1]
                temp_end = CNVscaffold[i, 2]
                # Ensure indices are within valid range
                frag_start = max(1, temp_star)
                frag_end = min(num_row, temp_end)
                if frag_start <= frag_end
                    depthFrag[i, j] = mean(AllMeanSub[frag_start:frag_end, j])
                else
                    depthFrag[i, j] = 0.0
                end
            end
        end
        
        # 3D table, 1D samples, 2D fragments, 3D copy number details
        copyNumberInfo = zeros(num_col, size(CNVscaffold, 1), length(gmmParam[1].means) + 1)
        for i in 1:num_col
            for j in 1:size(CNVscaffold, 1)
                try
                    # Use custom function to identify copy number
                    copyNumberInfo[i, j, :] = identify_copy_number(depthFrag[j, i], gmmParam[i], thresholdBFactor)
                    # Handle 0 copy number case
                    if depthFrag[j, i] < (gmmParam[i].means[1] - 2 * gmmParam[i].stds[1])
                        copyNumberInfo[i, j, 1] = 0
                    end
                catch e
                    @warn "Copy number identification failed for sample $i, fragment $j: " * string(e)
                    copyNumberInfo[i, j, 1] = -1  # Mark as undetermined
                end
            end
        end
        
        return copyNumberInfo
        
    catch e
        @error "Error occurred in copynumber_determination function: " * string(e)
        rethrow(e)
    end
end


"""
    genome_segment(inputDataDF, givenDistParam, breaks_threshold, flanking_length, minimum_length)

Identify genomic fragment breakpoints and construct CNV examination intervals

# Arguments:
- `inputDataDF`: Filtered and normalized data frame from function filtering_normalizing_seqDepth()
- `givenDistParam`: Distribution parameters for all gene fragments from function fitting_gmm()
- `breaks_threshold`: Breakpoint threshold, range [0, 1], e.g., 0.001
- `flanking_length`: Left and right flanking boundary length of examined fragment, should be consistent with GMM fitting fragment length, e.g., 10
- `minimum_length`: Minimum fragment length, e.g., 10

# Dependencies:
- rate_coincident: Calculate consistency probability between data and distribution parameters
- run_length_nonzero: Calculate non-zero run lengths

# Returns:
- `CNVscaffold`: CNV examination interval scaffold matrix

# Exceptions:
- Throws ArgumentError if input parameter types are incorrect
- Throws ArgumentError if parameter values exceed valid range
"""
function genome_segment(inputDataDF, givenDistParam, breaks_threshold, flanking_length, minimum_length)
    try
        # Input validation
        if !isa(inputDataDF, DataFrame)
            throw(ArgumentError("inputDataDF must be DataFrame type"))
        end
        
        if !isa(givenDistParam, Vector)
            throw(ArgumentError("givenDistParam must be Vector type"))
        end
        
        if !isa(breaks_threshold, Number) || breaks_threshold < 0 || breaks_threshold > 1
            throw(ArgumentError("breaks_threshold must be a number in [0, 1] range"))
        end
        
        if !isa(flanking_length, Integer) || flanking_length <= 0
            throw(ArgumentError("flanking_length must be a positive integer"))
        end
        
        if !isa(minimum_length, Integer) || minimum_length <= 0
            throw(ArgumentError("minimum_length must be a positive integer"))
        end
        
        # Check DataFrame has enough columns (at least 4: chr, binLowEdge, binUpEdge, depth data)
        if size(inputDataDF, 2) < 4
            throw(DimensionMismatch("inputDataDF requires at least 4 columns"))
        end
        
        # ***************Find all candidate breakpoints (CNV fragment boundaries)*******************
        subData = Matrix(inputDataDF[:, 4:end])
        row_num = size(subData, 1)
        sample_num = size(subData, 2)
        
        # Check parameter compatibility
        if flanking_length * 2 >= row_num
            throw(ArgumentError("flanking_length too large, incompatible with data length"))
        end
        
        if length(givenDistParam) != sample_num
            throw(DimensionMismatch("givenDistParam length must match sample count"))
        end
        
        LeftMeanSub = zeros(size(subData))
        RightMeanSub = zeros(size(subData))
        for i in (flanking_length + 1):(row_num - flanking_length)
            for j in 1:sample_num
                LeftMeanSub[i, j] = mean(subData[(i - flanking_length):(i - 1), j])
                RightMeanSub[i, j] = mean(subData[(i + 1):(i + flanking_length), j])
            end
        end
        
        temp_breaks = zeros(size(subData))
        for i in (flanking_length + 1):(row_num - flanking_length)
            for j in 1:sample_num
                temp = [LeftMeanSub[i, j]; RightMeanSub[i, j]]
                # Determine candidate breakpoints via probability from custom function
                if (rate_coincident(temp, givenDistParam[j]) < breaks_threshold) || 
                    (inputDataDF[i + flanking_length, 3] - inputDataDF[i - flanking_length, 3] > 100*10^3) # Large sequencing depth difference or physical distance too long, set as breakpoint
                    temp_breaks[i, j] = 1
                end
            end
        end
        
        temp = sum(temp_breaks, dims = 2) # Cumulative break count for each fragment across all samples
        
        # **************Consider run lengths between breakpoints********************
        temp_runs = run_length_nonzero(temp) # Merge breakpoints from all samples
        
        # Use intervals between adjacent break blocks as basic unit, intervals include boundaries
        # Calculate intervals
        if size(temp_runs, 1) <= 1
            @warn "Insufficient breakpoint runs detected"
            return zeros(0, 3)
        end
        
        temp_runs_region = zeros(size(temp_runs, 1) - 1, 3)
        for i in 1:(size(temp_runs, 1) - 1)
            temp_star = temp_runs[i, 3]
            temp_end = temp_runs[i + 1, 2]
            temp_len = temp_end - temp_star
            temp_runs_region[i, :] = [temp_star temp_end temp_len] 
        end
        
        # **Select intervals with fragment length >= minimum_length as final CNV examination intervals****
        temp_indices = findall(temp_runs_region[:, 3] .>= minimum_length)
        
        if isempty(temp_indices)
            @warn "No CNV intervals meeting minimum length requirement found"
            return zeros(0, 3)
        end
        
        # Filter out cross-chromosome fragments to form basic examination scaffold
        tempCNVscaffold = Int64.(temp_runs_region[temp_indices, :]) 
        temp_indices = findall(inputDataDF[tempCNVscaffold[:, 1], 1] .== inputDataDF[tempCNVscaffold[:, 2], 1])
        CNVscaffold = tempCNVscaffold[temp_indices, :]
        
        return CNVscaffold
        
    catch e
        @error "Error occurred in genome_segment function: " * string(e)
        rethrow(e)
    end
end


"""
    fitting_gmm(inputData, pre_mean, windowSize)

Fit data using Gaussian Mixture Model (GMM)

# Arguments:
- `inputData`: M×N matrix, N is sample count
- `pre_mean`: Predefined means for each component in Gaussian mixture model
- `windowSize`: Window size

# Returns:
- `EstimateParam`: Array containing GMM parameter estimation results for each sample

# Exceptions:
- Throws ArgumentError if input parameter types are incorrect
- Throws DimensionMismatch if inputData dimensions don't match windowSize
"""
function fitting_gmm(inputData, pre_mean, windowSize)
    try
        # Input validation
        if !isa(inputData, Matrix)
            throw(ArgumentError("inputData must be Matrix type"))
        end
        
        if !isa(pre_mean, Vector)
            throw(ArgumentError("pre_mean must be Vector type"))
        end
        
        if !isa(windowSize, Integer) || windowSize <= 0
            throw(ArgumentError("windowSize must be a positive integer"))
        end
        
        # Check dimension compatibility
        M, N = size(inputData)
        if windowSize > M
            throw(DimensionMismatch("windowSize cannot exceed inputData row count"))
        end
        
        # Initialize window mean matrix
        windowMean = zeros(M - windowSize + 1, N)
        
        # Calculate sliding window means
        for i in 1:(M - windowSize + 1)
            for j in 1:N
                # Ensure indices don't exceed bounds
                startIdx = max(1, i)
                endIdx = min(M, i + windowSize - 1)
                windowMean[i, j] = mean(inputData[startIdx:endIdx, j])
            end
        end
        
        # Get component count
        num_components = length(pre_mean)
        
        # If component count doesn't match predefined means count, give warning
        if num_components <= 0
            @warn "Predefined means count is zero"
            return []
        end
        
        # Parameter estimation
        EstimateParam = []
        for i in 1:N
            try
                # Use estimate_gmm function for parameter estimation
                gmm_result = estimate_gmm(windowMean[:, i], num_components; means = pre_mean)
                EstimateParam = [EstimateParam; gmm_result]
                print("Processing sample $i/$N\n")
            catch e
                @warn "Parameter estimation failed for sample $i: " * string(e)
            end
        end
        
        return EstimateParam
        
    catch e
        @error "Error occurred in fitting_gmm function: " * string(e)
        rethrow(e)
    end
end


"""
    filtering_normalizing_seqDepth(dataframe_raw_sequencingdepth, depth_filter1, depth_filter2)

Filter and normalize sequencing depth data

The main purpose of this function is to filter out fragments with abnormally low or high sequencing depth
due to various reasons, as these fragments have questionable data reliability, and normalize the remaining data.

# Arguments:
- `dataframe_raw_sequencingdepth`: DataFrame object containing sequencing depth information
    - Column 1: Chromosome ID (chr)
    - Column 2: Bin start position (binLowEdge)
    - Column 3: Bin end position (binUpEdge)
    - Column 4 and beyond: Coverage depth for each sample in this interval
- `depth_filter1`: Vector of length 2, representing min and max total depth thresholds before normalization
    - Element 1: Minimum depth threshold
    - Element 2: Maximum depth threshold
- `depth_filter2`: Vector of length 2, representing min and max total depth thresholds after normalization
    - Element 1: Minimum depth threshold
    - Element 2: Maximum depth threshold

# Returns:
- `normalSubCombineData`: Normalized and filtered sequencing depth data

# Exceptions:
- Throws ArgumentError if input parameter types are incorrect
- Throws DimensionMismatch if depth_filter length is not 2
"""
function filtering_normalizing_seqDepth(dataframe_raw_sequencingdepth, depth_filter1, depth_filter2)
    try
        # Input validation
        if !isa(dataframe_raw_sequencingdepth, DataFrame)
            throw(ArgumentError("dataframe_raw_sequencingdepth must be DataFrame type"))
        end
        
        if (!isa(depth_filter1, Vector) || length(depth_filter1) != 2) && (isa(depth_filter2, Vector) || length(depth_filter2) != 2)
            throw(DimensionMismatch("depth_filter must be a vector of length 2"))
        end
        
        if depth_filter1[1] >= depth_filter1[2]
            throw(ArgumentError("Minimum depth in depth_filter1 must be less than maximum depth"))
        end

        if depth_filter2[1] >= depth_filter2[2]
            throw(ArgumentError("Minimum depth in depth_filter2 must be less than maximum depth"))
        end
        
        # Check DataFrame has enough columns
        if size(dataframe_raw_sequencingdepth, 2) < 4
            throw(DimensionMismatch("dataframe_raw_sequencingdepth requires at least 4 columns"))
        end
        
        # Check total coverage depth
        tempDepth = sum(Matrix(dataframe_raw_sequencingdepth[:, 4:end]), dims = 2)
        
        # Extract only regions with sequencing depth in appropriate range
        # Extract fragments with total sequencing depth in specified range
        tempidx = findall((tempDepth[:] .> depth_filter1[1]) .& (tempDepth[:] .< depth_filter1[2]) .== 1)
        
        # If no rows meet conditions, return empty DataFrame
        if isempty(tempidx)
            @warn "No fragments meeting depth filter conditions"
            return DataFrame()
        end
        
        subCombineData = dataframe_raw_sequencingdepth[tempidx, :]
        
        # Normalize mean sequencing depth to 100
        tempMean = mean(Matrix(subCombineData[:, 4:end]), dims = 1)
        normalSubCombineData = hcat(subCombineData[:, 1:3], Float64.(subCombineData[:, 4:end]))
        
        # Normalize each row
        for i in 1:size(normalSubCombineData, 1)
            # Prevent division by zero error
            if all(tempMean .> 0)
                normalSubCombineData[i, 4:end] = (Vector{Float64}(normalSubCombineData[i, 4:end])' ./ tempMean) .* 100
            else
                @warn "Mean depth contains zero values, skipping normalization"
                break
            end
        end

        # Extract subset where min and max values of each row meet requirements and return
        miniDepth = minimum(Matrix(normalSubCombineData[:,4:end]), dims = 2)
        maxiDepth = maximum(Matrix(normalSubCombineData[:,4:end]), dims = 2)
        tempidx = findall((miniDepth[:] .> depth_filter2[1]) .& (maxiDepth[:] .< depth_filter2[2]))
        
        return normalSubCombineData[tempidx,:]
        
    catch e
        @error "Error occurred in filtering_normalizing_seqDepth function: " * string(e)
        rethrow(e)
    end
end
