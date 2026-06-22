% Distributions are created of the volumes for the smoothed model, TreeQSM and FuzzyQSM and their similarity is determined

clear variables
close all
clc

%% Inputs %%
    % Data locations
    smoothed_models_folder      =   'Data/Smoothed_Models';                         % Folder in which the smoothed models are located
    TreeQSM_folder_1            =   'Data/TreeQSMs/0175_beamdivergence';            % Folder in which the TreeQSM generated QSMs are stored
    TreeQSM_folder_2            =   'Data/TreeQSMs/030_beamdivergence';             % Folder in which the TreeQSM generated QSMs are stored
    FuzzyQSM_folder_1           =   'Data/FuzzyQSMs/0175_beamdivergence';           % Folder in which the FuzzyQSM generated QSMs are stored
    FuzzyQSM_folder_2           =   'Data/FuzzyQSMs/030_beamdivergence';            % Folder in which the FuzzyQSM generated QSMs are stored
    output_folder               =   'Data/Results/Volume_Distributions';            % Folder in which the results will be saved

    % Binning
    bin_diameter_width          =   0.01;                           % [m] Width of each bin in terms of diameter
        
    % Outputs
    num_decimal_digits          =   4;                              % [-]
    output_format               =   'Float';                        % [Float, Integer, Exponential]
    Print                       =   true;                           % [true, false] Whether or not the tables are displayed

%% Data %%
    % Beam divergences
    TreeQSM_folders         = {TreeQSM_folder_1, TreeQSM_folder_2};
    FuzzyQSM_folders        = {FuzzyQSM_folder_1, FuzzyQSM_folder_2};
    number_beam_divergences = length(TreeQSM_folders);

    beam_divergence_list = zeros(1, number_beam_divergences);

    for b = 1 : number_beam_divergences
        folder_name     = TreeQSM_folders{b};
        name_parts      = strsplit(folder_name, '/');

        Beam_Div_Part_fun   = @(part) contains(part, 'beamdivergence');
        beam_div_bool       = cellfun(Beam_Div_Part_fun, name_parts);
        beam_div_part       = name_parts{beam_div_bool};

        beam_divergence = strrep(beam_div_part, '_beamdivergence', '');
        beam_divergence = sprintf('%s.%s', beam_divergence(1), beam_divergence(2 : end));
        beam_divergence = str2double(beam_divergence);

        beam_divergence_list(b) = beam_divergence;
    end

    % Analysed trees
    smoothed_model_folders  = dir(smoothed_models_folder);
    smoothed_model_folders  = {smoothed_model_folders(:).name};

    dot_entries     = contains(smoothed_model_folders, '.');       % Remove the . and .. folders and the .xls if it has been created already
    tree_ID_cell    = smoothed_model_folders(~dot_entries);
    number_trees    = length(tree_ID_cell);

    % Determine the species of each of the trees
    Tree_ID_Split_fun   = @(tree_ID) strsplit(tree_ID, '_');
    tree_ID_parts_cell  = cellfun(Tree_ID_Split_fun, tree_ID_cell, 'UniformOutput', false);
    Tree_Species_fun    = @(tree_ID_parts) tree_ID_parts{1};
    tree_species_cell   = cellfun(Tree_Species_fun, tree_ID_parts_cell, 'UniformOutput', false);

    unique_tree_species = unique(tree_species_cell);
    number_species      = length(unique_tree_species);

    % Boolean for each species
    species_bool_cell = cell(1, number_species);
    
    for s = 1 : number_species
        tree_species            = unique_tree_species{s};
        Species_fun             = @(species) strcmp(species, tree_species);
        species_bool            = cellfun(Species_fun, tree_species_cell);
        species_bool_cell{s}    = species_bool;
    end

    % Including all trees
    bool_cell       = [{true(1, number_trees)}, species_bool_cell];
    bool_label_cell = ['All', unique_tree_species];
    number_bools    = length(bool_cell);

    % Volumes
    branch_segment_volume_cell      = cell(1, number_trees);
    branch_segment_min_radius_cell  = cell(1, number_trees);
    branch_segment_max_radius_cell  = cell(1, number_trees);
    TreeQSM_branch_cyl_volume_cell  = cell(number_beam_divergences, number_trees);
    TreeQSM_branch_cyl_radius_cell  = cell(number_beam_divergences, number_trees);
    FuzzyQSM_branch_cyl_volume_cell = cell(number_beam_divergences, number_trees);
    FuzzyQSM_branch_cyl_radius_cell = cell(number_beam_divergences, number_trees);

    for t = 1 : number_trees
        % Tree's ID
        tree_ID = tree_ID_cell{t};

        % Smoothed model segments
        smoothed_model_metrics_file_path    = sprintf('%s/%s/%s_Metrics.mat', smoothed_models_folder, tree_ID, tree_ID);
        Smoothed_Model_Metrics_File         = load(smoothed_model_metrics_file_path);
        
        Smoothed_Model_Circle_Geometry      = Smoothed_Model_Metrics_File.Smoothed_Model_Circle_Geometry;
        circle_radius_list                  = Smoothed_Model_Circle_Geometry.radius;
        
        Smoothed_Model_Metrics              = Smoothed_Model_Metrics_File.Smoothed_Model_Metrics;
        segment_branch_order_list           = Smoothed_Model_Metrics.segment.BranchOrder;
        segment_branch_bool                 = segment_branch_order_list ~= 0;

        branch_segment_volume_list          = Smoothed_Model_Metrics.segment.volume(segment_branch_bool);
        branch_segment_volume_cell{t}       = branch_segment_volume_list;

        branch_segment_circle_indices       = Smoothed_Model_Metrics.segment.circle_indices(segment_branch_bool);
        branch_segment_circle_radius_matrix = circle_radius_list(branch_segment_circle_indices);
        branch_segment_min_radius_list      = min(branch_segment_circle_radius_matrix, [], 2);
        branch_segment_min_radius_cell{t}   = double(branch_segment_min_radius_list);
        branch_segment_max_radius_list      = max(branch_segment_circle_radius_matrix, [], 2);
        branch_segment_max_radius_cell{t}   = double(branch_segment_max_radius_list);

        for b = 1 : number_beam_divergences
            % TreeQSM cylinder volumes
            TreeQSM_folder      = TreeQSM_folders{b};
            TreeQSM_file_path   = sprintf('%s/%s_Discrete_QSM.mat', TreeQSM_folder, tree_ID); 
            TreeQSM_File        = load(TreeQSM_file_path);
            QSM_Discrete        = TreeQSM_File.QSM_Discrete;

            TreeQSM_cyl_radius_list                 = QSM_Discrete.cylinder.radius;
            TreeQSM_cyl_length_list                 = QSM_Discrete.cylinder.length;
            TreeQSM_cyl_volume_list                 = pi * TreeQSM_cyl_radius_list.^2 .* TreeQSM_cyl_length_list;
            TreeQSM_cyl_BO_list                     = QSM_Discrete.cylinder.BranchOrder;
            fitted_stem_bool                        = TreeQSM_cyl_BO_list == 0;                         % Should be the same for FuzzyQSM
            TreeQSM_branch_cyl_volume_list          = TreeQSM_cyl_volume_list(~fitted_stem_bool);
            TreeQSM_branch_cyl_volume_cell{b, t}    = double(TreeQSM_branch_cyl_volume_list);
            TreeQSM_branch_cyl_radius_list          = TreeQSM_cyl_radius_list(~fitted_stem_bool);
            TreeQSM_branch_cyl_radius_cell{b, t}    = double(TreeQSM_branch_cyl_radius_list);

            % FuzzyQSM cylinder volumes
            FuzzyQSM_folder     = FuzzyQSM_folders{b};
            FuzzyQSM_file_path  = sprintf('%s/%s_Fuzzy_QSM.mat', FuzzyQSM_folder, tree_ID); 
            FuzzyQSM_File       = load(FuzzyQSM_file_path);
            QSM_Fuzzy           = FuzzyQSM_File.QSM_Fuzzy;

            FuzzyQSM_cyl_volume_list                = QSM_Fuzzy.cylinder.volume;
            FuzzyQSM_branch_cyl_volume_list         = FuzzyQSM_cyl_volume_list(~fitted_stem_bool);
            FuzzyQSM_branch_cyl_volume_cell{b, t}   = FuzzyQSM_branch_cyl_volume_list;
            FuzzyQSM_cyl_radius_list                = QSM_Fuzzy.cylinder.radius;
            FuzzyQSM_branch_cyl_radius_list         = FuzzyQSM_cyl_radius_list(~fitted_stem_bool);
            FuzzyQSM_branch_cyl_radius_cell{b, t}   = FuzzyQSM_branch_cyl_radius_list;
        end
    end

    % All smoothed model data
    branch_segment_volume_list      = vertcat(branch_segment_volume_cell{:});
    branch_segment_min_radius_list  = vertcat(branch_segment_min_radius_cell{:});
    branch_segment_max_radius_list  = vertcat(branch_segment_max_radius_cell{:});
    number_branch_segments          = length(branch_segment_max_radius_list);

    % Defining radius bins
    TreeQSM_branch_cyl_radius_list  = vertcat(TreeQSM_branch_cyl_radius_cell{:});
    FuzzyQSM_branch_cyl_radius_list = vertcat(FuzzyQSM_branch_cyl_radius_cell{:});

    max_radius              = max([branch_segment_max_radius_list; TreeQSM_branch_cyl_radius_list; FuzzyQSM_branch_cyl_radius_list]);
    bin_radius_edges_list   = 0 : bin_diameter_width/2 : max_radius + bin_diameter_width/2;
    number_bins             = length(bin_radius_edges_list) - 1;

%% Relative volume distributions for individual trees %%
    % Branch volumes are binned in terms of diameter and normalised by the true branch volume
    TreeQSM_RMSE_matrix         = zeros(number_beam_divergences, number_trees);
    FuzzyQSM_RMSE_matrix        = zeros(number_beam_divergences, number_trees);

    segment_rel_vol_distr_cell  = cell(1, number_trees);
    TreeQSM_rel_vol_distr_cell  = cell(number_beam_divergences, number_trees);
    FuzzyQSM_rel_vol_distr_cell = cell(number_beam_divergences, number_trees);

    for t = 1 : number_trees
        % Smoothed model data
        branch_segment_volume_list      = branch_segment_volume_cell{t};
        branch_segment_min_radius_list  = branch_segment_min_radius_cell{t};
        branch_segment_max_radius_list  = branch_segment_max_radius_cell{t};
        number_branch_segments          = length(branch_segment_max_radius_list);

        % Binning the segments
        bin_segment_volume_list = zeros(number_bins, 1);

        for i = 1 : number_bins
            % Bin extents
            bin_radius_min = bin_radius_edges_list(i);
            bin_radius_max = bin_radius_edges_list(i + 1);

            %--% Fraction of volume inside the bin %--%
            volume_fraction_list = zeros(number_branch_segments, 1);

            % The segment is fully contained within the bin
            fully_contained_bool                        = branch_segment_min_radius_list >= bin_radius_min & branch_segment_max_radius_list <= bin_radius_max;
            volume_fraction_list(fully_contained_bool)  = 1;               
    
            % The segment 'spills out' of the bin on both sides
            spilling_out_bool                           = branch_segment_min_radius_list < bin_radius_min & branch_segment_max_radius_list > bin_radius_max;
            volume_fraction_list(spilling_out_bool)     = (bin_radius_max - bin_radius_min) ./ (branch_segment_max_radius_list(spilling_out_bool) - branch_segment_min_radius_list(spilling_out_bool));
    
            % The segment spills out of the minimum side
            spilling_out_min_bool                       = branch_segment_min_radius_list < bin_radius_min & branch_segment_max_radius_list > bin_radius_min;
            spilling_out_min_bool                       = spilling_out_min_bool & ~spilling_out_bool;
            volume_fraction_list(spilling_out_min_bool) = (branch_segment_max_radius_list(spilling_out_min_bool) - bin_radius_min) ./ (branch_segment_max_radius_list(spilling_out_min_bool) - branch_segment_min_radius_list(spilling_out_min_bool));
    
            % The segment spills out of the maximum side
            spilling_out_max_bool                       = branch_segment_min_radius_list < bin_radius_max & branch_segment_max_radius_list > bin_radius_max;
            spilling_out_max_bool                       = spilling_out_max_bool & ~spilling_out_bool;
            volume_fraction_list(spilling_out_max_bool) = (bin_radius_max - branch_segment_min_radius_list(spilling_out_max_bool)) ./ (branch_segment_max_radius_list(spilling_out_max_bool) - branch_segment_min_radius_list(spilling_out_max_bool));
            
            %--% Total volume within the bin %--%
            bin_segment_volume_list(i) = sum(volume_fraction_list .* branch_segment_volume_list);                
        end
        
        segment_branch_volume           = sum(branch_segment_volume_list);
        bin_segment_rel_volume_list     = bin_segment_volume_list / segment_branch_volume;
        segment_rel_vol_distr_cell{t}   = bin_segment_rel_volume_list;

        for b = 1 : number_beam_divergences
            % This beam divergence's cylinder volumes and radii
            TreeQSM_branch_cyl_volume_list  = TreeQSM_branch_cyl_volume_cell{b, t};
            TreeQSM_branch_cyl_radius_list  = TreeQSM_branch_cyl_radius_cell{b, t};
            FuzzyQSM_branch_cyl_volume_list = FuzzyQSM_branch_cyl_volume_cell{b, t};
            FuzzyQSM_branch_cyl_radius_list = FuzzyQSM_branch_cyl_radius_cell{b, t};

            % Binning the QSM cylinder volumes
            [TreeQSM_bin_avg_volume_list, TreeQSM_bin_number_cyls_list, ~]      = Data_Binning(TreeQSM_branch_cyl_volume_list, TreeQSM_branch_cyl_radius_list, bin_radius_edges_list, number_bins);
            TreeQSM_bin_volume_list                                             = TreeQSM_bin_avg_volume_list .* TreeQSM_bin_number_cyls_list;
            TreeQSM_bin_rel_volume_list                                         = TreeQSM_bin_volume_list / segment_branch_volume;                  % Normalised by the true branch volume
            TreeQSM_bin_rel_volume_list(isnan(TreeQSM_bin_rel_volume_list))     = 0;
            TreeQSM_rel_vol_distr_cell{b, t}                                    = TreeQSM_bin_rel_volume_list;

            [FuzzyQSM_bin_avg_volume_list, FuzzyQSM_bin_number_cyls_list, ~]    = Data_Binning(FuzzyQSM_branch_cyl_volume_list, FuzzyQSM_branch_cyl_radius_list, bin_radius_edges_list, number_bins);
            FuzzyQSM_bin_volume_list                                            = FuzzyQSM_bin_avg_volume_list .* FuzzyQSM_bin_number_cyls_list;
            FuzzyQSM_bin_rel_volume_list                                        = FuzzyQSM_bin_volume_list / segment_branch_volume;                 % Normalised by the true branch volume
            FuzzyQSM_bin_rel_volume_list(isnan(FuzzyQSM_bin_rel_volume_list))   = 0;
            FuzzyQSM_rel_vol_distr_cell{b, t}                                   = FuzzyQSM_bin_rel_volume_list;

            % RMSE
            TreeQSM_SE                  = (bin_segment_rel_volume_list - TreeQSM_bin_rel_volume_list).^2;
            TreeQSM_RMSE                = sqrt(mean(TreeQSM_SE));
            TreeQSM_RMSE_matrix(b, t)   = TreeQSM_RMSE;

            FuzzyQSM_SE                 = (bin_segment_rel_volume_list - FuzzyQSM_bin_rel_volume_list).^2;
            FuzzyQSM_RMSE               = sqrt(mean(FuzzyQSM_SE));
            FuzzyQSM_RMSE_matrix(b, t)  = FuzzyQSM_RMSE; 
        end
    end

%% Relative volume distributions for all trees combined %%
    % True segment volumes
    total_segment_rel_vol_distr_matrix  = horzcat(segment_rel_vol_distr_cell{:});
    total_segment_rel_vol_distr_list    = sum(total_segment_rel_vol_distr_matrix, 2) / number_trees;

    % Per branch order
    total_TreeQSM_RMSE_list             = zeros(number_beam_divergences, 1);
    total_TreeQSM_rel_vol_distr_cell    = cell(number_beam_divergences, 1);
    total_FuzzyQSM_RMSE_list            = zeros(number_beam_divergences, 1);
    total_FuzzyQSM_rel_vol_distr_cell   = cell(number_beam_divergences, 1);

    for b = 1 : number_beam_divergences
        % TreeQSM
        total_TreeQSM_rel_vol_distr_matrix  = horzcat(TreeQSM_rel_vol_distr_cell{b, :});
        total_TreeQSM_rel_vol_distr_list    = sum(total_TreeQSM_rel_vol_distr_matrix, 2) / number_trees;
        total_TreeQSM_rel_vol_distr_cell{b} = total_TreeQSM_rel_vol_distr_list;

        total_TreeQSM_SE                    = (total_TreeQSM_rel_vol_distr_list - total_segment_rel_vol_distr_list).^2;
        total_TreeQSM_RMSE                  = sqrt(mean(total_TreeQSM_SE));
        total_TreeQSM_RMSE_list(b)          = total_TreeQSM_RMSE;

        % FuzzyQSM
        total_FuzzyQSM_rel_vol_distr_matrix     = horzcat(FuzzyQSM_rel_vol_distr_cell{b, :});
        total_FuzzyQSM_rel_vol_distr_list       = sum(total_FuzzyQSM_rel_vol_distr_matrix, 2) / number_trees;
        total_FuzzyQSM_rel_vol_distr_cell{b}    = total_FuzzyQSM_rel_vol_distr_list;
        total_FuzzyQSM_SE                       = (total_FuzzyQSM_rel_vol_distr_list - total_segment_rel_vol_distr_list).^2;
        total_FuzzyQSM_RMSE                     = sqrt(mean(total_FuzzyQSM_SE));
        total_FuzzyQSM_RMSE_list(b)             = total_FuzzyQSM_RMSE;
    end

    % Delta
    Delta_fun                               = @(total_rel_vol_distr_list) total_rel_vol_distr_list - total_segment_rel_vol_distr_list;
    delta_total_TreeQSM_rel_vol_distr_cell  = cellfun(Delta_fun, total_TreeQSM_rel_vol_distr_cell, 'UniformOutput', false);
    delta_total_FuzzyQSM_rel_vol_distr_cell = cellfun(Delta_fun, total_FuzzyQSM_rel_vol_distr_cell, 'UniformOutput', false);

    % Cumulative
    total_segment_cumul_rel_vol_distr_list  = cumsum(total_segment_rel_vol_distr_list);
    total_TreeQSM_cumul_rel_vol_distr_cell  = cellfun(@cumsum, total_TreeQSM_rel_vol_distr_cell, 'UniformOutput', false);
    total_FuzzyQSM_cumul_rel_vol_distr_cell = cellfun(@cumsum, total_FuzzyQSM_rel_vol_distr_cell, 'UniformOutput', false);

%% Plot %%
    % Method colours
    set_2_cmap      = cbrewer('qual', 'Set2', 3);
    TreeQSM_colour  = set_2_cmap(1, :);
    set_1_cmap      = cbrewer('qual', 'Set1', 3);
    FuzzyQSM_colour = set_1_cmap(1, :);

    % Plotted data sets
    data_labels         = {'relative volume', '\Delta relative volume', 'cumulative rel. volume'};
    data_units          = {'-', '-', '-'};
    number_data_sets    = length(data_units);

    segment_data_cell   = {total_segment_rel_vol_distr_list, [], total_segment_cumul_rel_vol_distr_list};
    TreeQSM_data_cell   = {total_TreeQSM_rel_vol_distr_cell, delta_total_TreeQSM_rel_vol_distr_cell, total_TreeQSM_cumul_rel_vol_distr_cell};
    FuzzyQSM_data_cell  = {total_FuzzyQSM_rel_vol_distr_cell, delta_total_FuzzyQSM_rel_vol_distr_cell, total_FuzzyQSM_cumul_rel_vol_distr_cell};

    % Plots
    for d = 1 : number_data_sets
        % The data set
        data_label          = data_labels{d};
        data_unit           = data_units{d};
        segment_data_list   = segment_data_cell{d};

        for b = 1 : number_beam_divergences
            % This beam divergence's data
            beam_divergence = beam_divergence_list(b);
    
            total_TreeQSM_rel_vol_distr_list    = TreeQSM_data_cell{d}{b};
            total_FuzzyQSM_rel_vol_distr_list   = FuzzyQSM_data_cell{d}{b};

            % The figure
            figure_number   = d*10 + b;
            Fig             = figure(figure_number);
            figure_name     = sprintf('%s distributions BD%.3g', data_label, beam_divergence);
            figure_name     = strrep(figure_name, '\', '');
            figure_name     = strrep(figure_name, '.', '');
        
            set(Fig, 'name', figure_name, 'NumberTitle', 'off');
        
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
        
            hold on
            grid on
        
            % Title
            sgtitle(sprintf('%s = %.3f mrad', '\lambda', beam_divergence), 'FontSize', 25);
                
            % True values
            if ~contains(data_label, '\Delta')
                bin_diameter_centre_list = (bin_radius_edges_list(1 : number_bins) + bin_radius_edges_list(2 : number_bins + 1));       % The factor 2 cancels out
                plot(bin_diameter_centre_list, segment_data_list, 'LineWidth', 3, 'Color', 'k', 'DisplayName', 'True');
            end
            
            % TreeQSM values
            plot(bin_diameter_centre_list, total_TreeQSM_rel_vol_distr_list, 'LineWidth', 3, 'Color', TreeQSM_colour, 'DisplayName', 'TreeQSM');
        
            % FuzzyQSM values
            plot(bin_diameter_centre_list, total_FuzzyQSM_rel_vol_distr_list, 'LineWidth', 3, 'Color', FuzzyQSM_colour, 'DisplayName', 'FuzzyQSM');
        
            % Axes
            xlim([0, 2*max_radius]);
            xlabel('diameter [m]');
            ylabel(sprintf('%s [%s]', data_label, data_unit));

            if strcmp(data_label, 'relative volume')
                ylim([0, 0.25]);
            elseif strcmp(data_label, '\Delta relative volume')
                ylim([-0.12, 0.12]);
            elseif strcmp(data_label, 'cumulative rel. volume')
                ylim([0, 1]);
            end
        
            % Formatting
            legend('show', 'location', 'eastoutside');
        
            set(gca, 'FontSize', 25);
            set(gca, 'LineWidth', 3);
        
            hold off
        
            % Saving the figure     
            export_fig(figure_number, sprintf('%s.png', figure_name));
            movefile([figure_name, '*'], output_folder); 
        end
    end

%% Table %%
    % Method indices for consistency
    TreeQSM_ind     = 1;
    FuzzyQSM_ind    = 2;
    number_methods  = FuzzyQSM_ind;

    % Average RMSE for each boolean
    mean_RMSE_matrix = zeros(number_methods * number_beam_divergences, number_bools);   

    for i = 1 : number_bools
        % This boolean's data
        bool_list       = bool_cell{i};
        
        mean_RMSE_list    = zeros(number_methods, number_beam_divergences);
        row_labels      = cell(number_methods, number_beam_divergences);

        for m = 1 : number_methods
            for b = 1 : number_beam_divergences
                beam_divergence = beam_divergence_list(b);

                if m == TreeQSM_ind
                    RMSE            = mean(TreeQSM_RMSE_matrix(b, bool_list));
                    method_label    = 'TreeQSM';
                elseif m == FuzzyQSM_ind
                    RMSE            = mean(FuzzyQSM_RMSE_matrix(b, bool_list));
                    method_label    = 'FuzzyQSM';
                end

                mean_RMSE_list(m, b)    = RMSE;
                row_labels{m, b}        = sprintf('%s, lambda = %.3f mrad', method_label, beam_divergence);
            end
        end

        mean_RMSE_matrix(:, i)  = mean_RMSE_list(:);
        row_labels              = reshape(row_labels, [number_methods * number_beam_divergences, 1]);
    end

    RMSE_table_file_name  = 'Relative_Volume_Distributions_RMSE.xls';
    Table_Formatter(mean_RMSE_matrix, num_decimal_digits, output_format, row_labels, bool_label_cell, RMSE_table_file_name, Print);

    movefile(RMSE_table_file_name, output_folder);
