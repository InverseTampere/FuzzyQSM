% After the metrics of the smoothed models of the Wytham trees have been determined, they are analysed here

clear variables
close all
beep off
clc

%% Inputs %%
    % Data locations
    smoothed_models_folder  = "C:\Users\svvive\OneDrive - TUNI.fi\Remote_Sensing_Forests\FuzzyQSM\Data\Smoothed_Models";
    original_QSMs_folder    = "C:\Users\svvive\OneDrive - TUNI.fi\Remote_Sensing_Forests\FuzzyQSM\Data\Original_QSMs";

    % Outputs
    num_decimal_digits          =   4;                      % [-]
    output_format               =   'Float';                % [Float, Integer, Exponential]
    Print                       =   true;                  % [true, false] Whether or not the tables are displayed

%% Analysed trees %%
    % Tree ID folders which contain the data
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

%% Tree metrics %%
    % Tree metrics structure
    Original_Tree_Metrics_cell = cell(1, number_trees);
    Smoothed_Tree_Metrics_cell = cell(1, number_trees);

    for t = 1 : number_trees
        % This tree's metrics
        tree_ID                     = tree_ID_cell{t};
        metrics_file_name           = sprintf("%s%s%s%s%s_Metrics.mat", smoothed_models_folder, '\', tree_ID, '\', tree_ID);
        Smoothed_Model_Metrics_File = load(metrics_file_name);

        %--% Original metrics %--%
        % Original QSM's
        Original_Tree_Metrics       = Smoothed_Model_Metrics_File.Original_QSM.treedata;
        Original_Branch_Metrics     = Smoothed_Model_Metrics_File.Original_QSM.branch; 
        Original_Cylinder_Metrics   = Smoothed_Model_Metrics_File.Original_QSM.cylinder; 

        % Conversion of the volumes from dm^3 to m^3
        Original_Tree_Metrics.TotalVolume   = Original_Tree_Metrics.TotalVolume / 1e3;                        
        Original_Tree_Metrics.BranchVolume  = Original_Tree_Metrics.BranchVolume / 1e3;                        
        Original_Tree_Metrics.TrunkVolume   = Original_Tree_Metrics.TrunkVolume / 1e3;                          

        % Computing the tree's total branch length
        branch_length_list                  = Original_Branch_Metrics.length;
        branch_order_list                   = Original_Branch_Metrics.order;
        branch_stem_bool                    = branch_order_list == 0;
        original_branch_length              = sum(branch_length_list(~branch_stem_bool));
        Original_Tree_Metrics.BranchLength  = original_branch_length;

        % Relative centre of gravity height
        cyl_radius_list     = Original_Cylinder_Metrics.radius;
        cyl_length_list     = Original_Cylinder_Metrics.length;
        cyl_volume_list     = pi*cyl_radius_list.^2 .* cyl_length_list;

        cyl_start_matrix    = Original_Cylinder_Metrics.start;
        cyl_axis_matrix     = Original_Cylinder_Metrics.axis;
        cyl_centre_matrix   = cyl_start_matrix + cyl_length_list/2 .* cyl_axis_matrix;
        cyl_height_list     = cyl_centre_matrix(:, 3);
        rel_cyl_height_list = cyl_height_list / max(cyl_height_list);

        original_rel_cog_height             = sum(rel_cyl_height_list .* cyl_volume_list) / sum(cyl_volume_list);
        Original_Tree_Metrics.RelCoGHeight  = original_rel_cog_height;

        cyl_branch_order_list   = Original_Cylinder_Metrics.BranchOrder;
        cyl_stem_bool           = cyl_branch_order_list == 0;

        original_rel_branch_cog_height              = sum(rel_cyl_height_list(~cyl_stem_bool) .* cyl_volume_list(~cyl_stem_bool)) / sum(cyl_volume_list(~cyl_stem_bool));
        Original_Tree_Metrics.RelBranchCoGHeight    = original_rel_branch_cog_height;

        % Adding to the cell array
        Original_Tree_Metrics_cell{t} = Original_Tree_Metrics;

        %--% Smoothed metrics %--%
        % Smoothed model's metrics
        Smoothed_Tree_Metrics       = Smoothed_Model_Metrics_File.Smoothed_Model_Metrics.treedata;
        Smoothed_Branch_Metrics     = Smoothed_Model_Metrics_File.Smoothed_Model_Metrics.branch;        
        Smoothed_Model_Segments     = Smoothed_Model_Metrics_File.Smoothed_Model_Metrics.segment;
        Smoothed_Model_Circles      = Smoothed_Model_Metrics_File.Smoothed_Model_Circle_Geometry;
            
        % Total branch length
        smoothed_branch_length_list         = Smoothed_Branch_Metrics.length;
        smoothed_branch_order_list          = Smoothed_Branch_Metrics.order;
        smoothed_stem_bool                  = smoothed_branch_order_list == 0;
        smoothed_branch_length              = sum(smoothed_branch_length_list(~smoothed_stem_bool));
        Smoothed_Tree_Metrics.BranchLength  = smoothed_branch_length;

        % Relative centre of gravity height
        segment_volume_list             = Smoothed_Model_Metrics_File.Smoothed_Model_Metrics.segment.volume;
        segment_circle_indices_matrix   = Smoothed_Model_Metrics_File.Smoothed_Model_Metrics.segment.circle_indices;

        circle_centre_matrix            = Smoothed_Model_Metrics_File.Smoothed_Model_Circle_Geometry.centre;
        circle_height_list              = circle_centre_matrix(:, 3);
        segment_circle_height_matrix    = circle_height_list(segment_circle_indices_matrix);
        segment_height_list             = mean(segment_circle_height_matrix, 2);
        segment_rel_height_list         = segment_height_list / max(segment_height_list);

        smoothed_rel_cog_height             = sum(segment_rel_height_list .* segment_volume_list) / sum(segment_volume_list);
        Smoothed_Tree_Metrics.RelCoGHeight  = smoothed_rel_cog_height;

        segment_branch_order_list   = Smoothed_Model_Metrics_File.Smoothed_Model_Metrics.segment.BranchOrder;
        segment_stem_bool           = segment_branch_order_list == 0;

        smoothed_rel_branch_cog_height              = sum(segment_rel_height_list(~segment_stem_bool) .* segment_volume_list(~segment_stem_bool)) / sum(segment_volume_list(~segment_stem_bool));
        Smoothed_Tree_Metrics.RelBranchCoGHeight    = smoothed_rel_branch_cog_height;

        % The structure is added to the cell array
        Smoothed_Tree_Metrics_cell{t} = Smoothed_Tree_Metrics;
    end

    % Structure containing the metrics
    Original_Tree_Metrics                                   = vertcat(Original_Tree_Metrics_cell{:});
    [original_tree_metrics_matrix, Original_Field_Indices]  = Structure_Data_Concatenation(Original_Tree_Metrics);

    Smoothed_Tree_Metrics                                   = vertcat(Smoothed_Tree_Metrics_cell{:});
    [smoothed_tree_metrics_matrix, Smoothed_Field_Indices]  = Structure_Data_Concatenation(Smoothed_Tree_Metrics);

    Tree_Metrics        = struct('Smoothed', struct('matrix', smoothed_tree_metrics_matrix, 'indices', Smoothed_Field_Indices, 'tree_metrics', {fieldnames(Smoothed_Field_Indices)}, 'folder', smoothed_models_folder), ...
                                 'Original', struct('matrix', original_tree_metrics_matrix, 'indices', Original_Field_Indices, 'tree_metrics', {fieldnames(Original_Field_Indices)}, 'folder', original_QSMs_folder));
    metric_sets         = fieldnames(Tree_Metrics);
    number_metric_sets  = length(metric_sets);

%% Tables %%
    % The mean, maximum and minimum are tabulated for each metric and for each species as well as the data as a whole
    bool_cell       = [{true(1, number_trees)}, species_bool_cell];
    bool_label_cell = ['All', unique_tree_species];
    number_bools    = length(bool_cell);

    % For each of the data sets they are generated separately
    for i = 1 : number_metric_sets
        % This set's data
        metric_set          = metric_sets{i};
        tree_metrics_matrix = Tree_Metrics.(metric_set).matrix;
        tree_metrics        = Tree_Metrics.(metric_set).tree_metrics;
        folder_name         = Tree_Metrics.(metric_set).folder;

        number_tree_metrics         = length(tree_metrics);
        tree_metrics_mean_matrix    = zeros(number_bools, number_tree_metrics);
        tree_metrics_max_matrix     = zeros(number_bools, number_tree_metrics);
        tree_metrics_min_matrix     = zeros(number_bools, number_tree_metrics);
    
        for b = 1 : number_bools
            % This boolean's data
            bool_list               = bool_cell{b};
            tree_metrics_matrix_b   = tree_metrics_matrix(bool_list, :);
    
            % Mean, max and min
            tree_metrics_mean_list          = mean(tree_metrics_matrix_b, 1);
            tree_metrics_mean_matrix(b, :)  = tree_metrics_mean_list;
            tree_metrics_max_list           = max(tree_metrics_matrix_b, [], 1);
            tree_metrics_max_matrix(b, :)   = tree_metrics_max_list;
            tree_metrics_min_list           = min(tree_metrics_matrix_b, [], 1);
            tree_metrics_min_matrix(b, :)   = tree_metrics_min_list;
        end
    
        % Tables
        mean_table_file_name = sprintf('%s_Tree_Metrics_Mean.xls', metric_set);
        Table_Formatter(tree_metrics_mean_matrix, num_decimal_digits, output_format, bool_label_cell, tree_metrics, mean_table_file_name, Print);
    
        max_table_file_name = sprintf('%s_Tree_Metrics_Max.xls', metric_set);
        Table_Formatter(tree_metrics_max_matrix, num_decimal_digits, output_format, bool_label_cell, tree_metrics, max_table_file_name, Print);
    
        min_table_file_name = sprintf('%s_Tree_Metrics_Min.xls', metric_set);
        Table_Formatter(tree_metrics_min_matrix, num_decimal_digits, output_format, bool_label_cell, tree_metrics, min_table_file_name, Print);
    
        % The tables are moved to the smoothed models folder
        movefile(sprintf('%s_Tree_Metrics*.xls', metric_set), folder_name);
    end

%% Plots %%
    % Species colour map
    species_cmap = cbrewer('qual', 'Set1', number_species);
    species_cmap = max(species_cmap, 0);
    species_cmap = min(species_cmap, 1);

    % Metrics are plotted versus the DBH
    DBH_field       = 'DBHqsm';
    figure_number   = 0;

    for i = 1 : number_metric_sets
        % This set's data
        metric_set          = metric_sets{i};
        tree_metrics_matrix = Tree_Metrics.(metric_set).matrix;
        Metrics_Indices     = Tree_Metrics.(metric_set).indices;
        tree_metrics        = Tree_Metrics.(metric_set).tree_metrics;
        folder_name         = Tree_Metrics.(metric_set).folder;

        number_tree_metrics = length(tree_metrics);

        DBH_index   = Metrics_Indices.(DBH_field);
        DBH_list    = tree_metrics_matrix(:, DBH_index);    
    
        for t = 1 : number_tree_metrics
            % This metric's data
            metric_field    = tree_metrics{t};
            metric_index    = Metrics_Indices.(metric_field);
            metric_list     = tree_metrics_matrix(:, metric_index);
    
            if contains(metric_field, 'DBH') 
                continue
            end
    
            % A neater label
            metric_label = regexprep(metric_field, '([A-Z])', ' ${lower($1)}');     % Inserts a space before each capital letter and makes them lower case
    
            if strcmp(metric_label(1), ' ')                                         % This may result in a space being placed at the start, which is removed
                metric_label(1) = [];
            end
    
            metric_label = strrep(metric_label, 'trunk', 'stem');                   % Trunk is replaced by stem
    
            % Unit
            if contains(metric_field, 'Volume')
                metric_unit     = 'm^3';
            elseif strcmp(metric_field, 'RelCoGHeight')
                metric_unit     = '-';
                metric_label    = 'rel. height c.o.g.';
            elseif strcmp(metric_field, 'RelBranchCoGHeight')
                metric_unit     = '-';
                metric_label    = 'rel. height branch c.o.g.';
            elseif contains(metric_field, 'Height') || contains(metric_field, 'Length')
                metric_unit     = 'm';
            end
    
            %--% Figure %--%
            figure_number   = figure_number + 1;
            Fig             = figure(figure_number);
            figure_name     = sprintf('%s %s vs. DBH', metric_set, metric_label);
    
            set(Fig, 'name', figure_name, 'NumberTitle', 'off');
    
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.85 0.5])
            set(gcf, 'color', [1, 1, 1])    
    
            hold on
            grid on
    
            % Each species gets its own colour
            for s = 1 : number_species
                % This species' data
                species         = unique_tree_species{s};
                species_colour  = species_cmap(s, :);
                species_bool    = species_bool_cell{s};
    
                % Scatter plot
                scatter(DBH_list(species_bool), metric_list(species_bool), 50, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', species_colour, 'DisplayName', species);
            end
    
            % Axes
            xlabel('DBH [m]');
            ylabel(sprintf('%s [%s]', metric_label, metric_unit));
    
            % Legend
            legend('show', 'location', 'eastoutside');
    
            set(gca, 'FontSize', 20);
            set(gca, 'LineWidth', 2);
    
            hold off
    
            % Saving the figure        
            export_fig(figure_number, sprintf('%s.fig', figure_name));
            export_fig(figure_number, sprintf('%s.png', figure_name));
    
            movefile([figure_name, '*'], folder_name); 
        end
    end



    