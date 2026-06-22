% This script tests fuzzy fitting only those cylinders with relative uncertainty greater than a given threshold
% Relative uncertainty is defined as sqrt(sigma_rad * sigma_prop) / cyl_radius
% Note that the FuzzyQSM and TreeQSM QSMs need to both have already been fitted

clear variables
close all
beep off
clc

%% Inputs %%
    %--% Relative uncertainty %--%
    max_threshold               =   1.0;                            % [-] Highest tested threshold 
    number_steps                =   011;                            % [-] Number of steps that are tested between 0 and the maximum

    %--% Data %--%
    smoothed_models_folder      =   'Data/Smoothed_Models';                                     % Folder in which the smoothed models are located
    TreeQSM_folder              =   'Data/TreeQSMs/030_beamdivergence';                         % Folder in which the TreeQSM generated QSMs are located
    FuzzyQSM_folder             =   'Data/FuzzyQSMs/030_beamdivergence';                        % Folder in which the TreeQSM generated QSMs are located
    output_folder               =   'Data/Results/030_beamdivergence/Relative_Uncertainty';     % Folder in which the results will be saved

    %--% Outputs %--%
    QSM_Print                   =   false;                  % [true, false] Shows printed statements regarding intermediate results
    QSM_Plot                    =   false;                  % [true, false] Shows plots of intermediate results
    Print                       =   true;                  % [true, false] Prints final results
    Plot                        =   true;                  % [true, false] Plots final results

    num_decimal_digits          =   4;                      % [-]
    output_format               =  'Float';                 % [Float, Integer, Exponential]
    
%% Analysed trees %%
    % Tree IDs
    TreeQSM_files   = dir(TreeQSM_folder);
    TreeQSM_files   = {TreeQSM_files(:).name};

    dot_folders     = startsWith(TreeQSM_files, '.');       % Remove the . and .. folders
    TreeQSM_files   = TreeQSM_files(~dot_folders);

    Tree_ID_fun     = @(TreeQSM_file) strrep(TreeQSM_file, '_Discrete_QSM.mat', '');
    tree_ID_cell    = cellfun(Tree_ID_fun, TreeQSM_files, 'UniformOutput', false);

    number_trees = length(tree_ID_cell);

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

%% Relative uncertainty effects %%
    % Relative uncertainty thresholds
    rel_uncertainty_threshold_steps     = linspace(0, max_threshold, number_steps);
    rel_uncertainty_threshold_labels    = string(rel_uncertainty_threshold_steps);

    % Cycling through the trees
    Results_Tables_cell             = cell(number_trees, number_steps);
    fuzzy_cyl_percentage_matrix     = zeros(number_trees, number_steps);
    fuzzy_points_percentage_matrix  = zeros(number_trees, number_steps);

    for t = 1 : number_trees
        % This tree's ID
        tree_ID = tree_ID_cell{t};

        % The TreeQSM fitted QSM
        TreeQSM_file            = TreeQSM_files{t};
        TreeQSM_file_name       = sprintf('%s/%s', TreeQSM_folder, TreeQSM_file);
        Tree_QSM_File           = load(TreeQSM_file_name);
        Discrete_QSM            = Tree_QSM_File.QSM_Discrete;
        TreeQSM_Inputs          = Tree_QSM_File.TreeQSM_Inputs;
        
        Discrete_Cylinders      = Discrete_QSM.cylinder;

        % The fuzzy QSM
        FuzzyQSM_file_name  = sprintf('%s/%s_Fuzzy_QSM.mat', FuzzyQSM_folder, tree_ID);
        Fuzzy_QSM_File      = load(FuzzyQSM_file_name);
        Fuzzy_QSM           = Fuzzy_QSM_File.QSM_Fuzzy;

        Fuzzy_Cylinders     = Fuzzy_QSM.cylinder;

        % Ensure the fields are compatible
        number_cylinders = length(Discrete_Cylinders.radius);

        Discrete_Cylinders.centre       = Discrete_Cylinders.start + Discrete_Cylinders.length/2 .* Discrete_Cylinders.axis;
        Discrete_Cylinders.volume       = pi * Discrete_Cylinders.radius.^2 .* Discrete_Cylinders.length;
        Discrete_Cylinders.fuzzy_bool   = false(number_cylinders, 1);
        Discrete_Cylinders.AxisChanged  = false(number_cylinders, 1);

        Fuzzy_Cylinders.extension       = Discrete_Cylinders.extension;
        Fuzzy_Cylinders.added           = Discrete_Cylinders.added;
        Fuzzy_Cylinders.mad             = Discrete_Cylinders.mad;

        % The point cloud distributions
        point_cloud_file_name               = sprintf('%s/%s_Point_Cloud_Distributions.mat', FuzzyQSM_folder, tree_ID);
        Point_Cloud_File                    = load(point_cloud_file_name);
        Cyl_Point_Cloud_Distributions_cell  = Point_Cloud_File.Cyl_Point_Cloud_Distributions_cell;
    
        empty_bool = cellfun(@isempty, Cyl_Point_Cloud_Distributions_cell);

        Number_Points_fun               = @(Point_Cloud_Distributions) length(Point_Cloud_Distributions.sigma_radial_cell);
        number_points_list              = zeros(1, number_cylinders);
        number_points_list(~empty_bool) = cellfun(Number_Points_fun, Cyl_Point_Cloud_Distributions_cell(~empty_bool));

        Point_Cloud_fun                         = @(Point_Cloud_Distributions) vertcat(Point_Cloud_Distributions.distribution_mu_cell{:});
        cylinder_point_cloud_cell               = cell(1, number_cylinders); 
        cylinder_point_cloud_cell(~empty_bool)  = cellfun(Point_Cloud_fun, Cyl_Point_Cloud_Distributions_cell(~empty_bool), 'UniformOutput', false);

        % Propagation and radial uncertainty
        Sigma_Radial_fun                = @(Cyl_Point_Cloud_Distributions) vertcat(Cyl_Point_Cloud_Distributions.sigma_radial_cell{:});
        sigma_radial_cell               = cell(1, number_cylinders);
        sigma_radial_cell(~empty_bool)  = cellfun(Sigma_Radial_fun, Cyl_Point_Cloud_Distributions_cell(~empty_bool), 'UniformOutput', false);
    
        Sigma_Prop_fun                  = @(Cyl_Point_Cloud_Distributions) vertcat(Cyl_Point_Cloud_Distributions.sigma_prop_cell{:});
        sigma_prop_cell                 = cell(1, number_cylinders);
        sigma_prop_cell(~empty_bool)    = cellfun(Sigma_Prop_fun, Cyl_Point_Cloud_Distributions_cell(~empty_bool), 'UniformOutput', false);
    
        % Geometric average uncertainty
        Geom_Average_Uncertainty_fun            = @(sigma_radial_list, sigma_prop_list) sqrt(sigma_radial_list .* sigma_prop_list);
        geom_avg_uncertainty_cell               = cellfun(Geom_Average_Uncertainty_fun, sigma_radial_cell, sigma_prop_cell, 'UniformOutput', false);
        geom_avg_uncertainty_list               = NaN(1, number_cylinders);
        geom_avg_uncertainty_list(~empty_bool)  = cellfun(@mean, geom_avg_uncertainty_cell(~empty_bool));

        % Relative uncertainty magnitude
        discr_cyl_radius_list       = Discrete_Cylinders.radius;
        rel_uncertainty_magn_list   = geom_avg_uncertainty_list ./ discr_cyl_radius_list';

        % Smoothed model metrics
        metrics_file_name               = sprintf('%s/%s/%s_Metrics.mat', smoothed_models_folder, tree_ID, tree_ID);
        Smoothed_Model_Metrics_File     = load(metrics_file_name);
        Smoothed_Model_Metrics          = Smoothed_Model_Metrics_File.Smoothed_Model_Metrics;
        Original_QSM                    = Smoothed_Model_Metrics_File.Original_QSM;
        Smoothed_Model_Circle_Geometry  = Smoothed_Model_Metrics_File.Smoothed_Model_Circle_Geometry;

        Original_QSM.treedata.TreeHeight = Smoothed_Model_Metrics.tree.TreeHeight;

        % Uncertainty threshold results
        for i = 1 : number_steps
            % Cylinders which exceed the threshold
            rel_uncertainty_threshold   = rel_uncertainty_threshold_steps(i);
            fuzzy_cyl_bool              = rel_uncertainty_magn_list > rel_uncertainty_threshold;

            % Percentage of fuzzy fitted cylinders
            fuzzy_cyl_percentage                = sum(fuzzy_cyl_bool) / number_cylinders * 100;
            fuzzy_cyl_percentage_matrix(t, i)   = fuzzy_cyl_percentage;
    
            % Percentage of points used for fuzzy fits
            fuzzy_points_percentage                 = sum(number_points_list(fuzzy_cyl_bool)) / sum(number_points_list) * 100;
            fuzzy_points_percentage_matrix(t, i)    = fuzzy_points_percentage;

            % Creating a mixed QSM
            cylinder_fields = fieldnames(Fuzzy_Cylinders);
            number_fields   = length(cylinder_fields);
            Mixed_Cylinders = struct();
    
            for f = 1 : number_fields
                % This field's data
                field               = cylinder_fields{f};
                fuzzy_field_data    = Fuzzy_Cylinders.(field);
                discr_field_data    = Discrete_Cylinders.(field);
    
                % Creating the mixed data
                num_dim             = size(fuzzy_field_data, 2);
                mixed_field_data    = zeros(number_cylinders, num_dim);
    
                mixed_field_data(fuzzy_cyl_bool, :)     = fuzzy_field_data(fuzzy_cyl_bool, :);
                mixed_field_data(~fuzzy_cyl_bool, :)    = discr_field_data(~fuzzy_cyl_bool, :);
    
                Mixed_Cylinders.(field) = mixed_field_data;
            end

            % Ensuring the branch cylinder axes are reasonable and point down the branch
            Mixed_Cylinders = Branch_Axis_Correction(Mixed_Cylinders, Discrete_Cylinders);
    
            % Correct the radii
            Mixed_Cylinders = Branch_Radius_Correction(Mixed_Cylinders, TreeQSM_Inputs, cylinder_point_cloud_cell);
    
            % Ensuring the cylinder interfaces are reasonable
            Mixed_Cylinders = Branch_Interface_Fitting(Mixed_Cylinders);

            % Branch metrics
            [Mixed_Branch_Metrics, Mixed_Branch_Order_Metrics] = QSM_Branch_Metrics(Mixed_Cylinders, Discrete_QSM);
    
            % Tree metrics
            Mixed_Tree_Metrics = QSM_Tree_Metrics(Mixed_Cylinders, Mixed_Branch_Order_Metrics);
    
            % Full QSM
            Mixed_QSM = struct('cylinder', Mixed_Cylinders, 'treedata', Mixed_Tree_Metrics, 'branch', Mixed_Branch_Metrics, 'branch_order', Mixed_Branch_Order_Metrics);

            % Results
            Output_Decisions            = struct('Print', QSM_Print, 'Plot', QSM_Plot, 'Compute_Radius_Errors', false);
            tree_ID_threshold           = sprintf('%s_%i', tree_ID, i);
            Results_Tables              = FuzzyQSM_Result_Evaluation(tree_ID_threshold, Original_QSM, Discrete_QSM, Mixed_QSM, Smoothed_Model_Metrics, Smoothed_Model_Circle_Geometry, Cyl_Point_Cloud_Distributions_cell, [], [], [], [], [], Output_Decisions);
            Results_Tables_cell{t, i}   = Results_Tables;
        end
    end

    % Save the results
    results_file_name = 'Rel_Uncertainty_Results.mat';
    save(results_file_name, 'Results_Tables_cell', 'fuzzy_cyl_percentage_matrix', 'fuzzy_points_percentage_matrix', 'rel_uncertainty_threshold_steps');
    movefile(results_file_name, output_folder);

%% Diagnostic metrics %%
    figure_number = 0;

    for b = 1 : number_bools
        % This boolean's data
        bool_list   = bool_cell{b};
        bool_label  = bool_label_cell{b};

        fuzzy_cyl_percentage_matrix_b       = fuzzy_cyl_percentage_matrix(bool_list, :);
        fuzzy_points_percentage_matrix_b    = fuzzy_points_percentage_matrix(bool_list, :);

        %--% Fuzzy cylinder percentage table %--%
        mean_fuzzy_cyl_perc_list    = mean(fuzzy_cyl_percentage_matrix_b, 1);
        std_fuzzy_cyl_perc_list     = std(fuzzy_cyl_percentage_matrix_b, [], 1);
        max_fuzzy_cyl_perc_list     = max(fuzzy_cyl_percentage_matrix_b, [], 1);
        min_fuzzy_cyl_perc_list     = min(fuzzy_cyl_percentage_matrix_b, [], 1);
    
        fuzzy_cyl_perc_table_matrix = [mean_fuzzy_cyl_perc_list; std_fuzzy_cyl_perc_list; max_fuzzy_cyl_perc_list; min_fuzzy_cyl_perc_list];
        row_labels                  = {'mean', 'std', 'max', 'min'};
    
        table_file_name = sprintf('Rel_Uncertainty_Fuzzy_Cyl_Percentage_%s.xls', bool_label);

        fprintf('Fuzzy cylinder percentage for %s trees \n', bool_label);
        Table_Formatter(fuzzy_cyl_perc_table_matrix, num_decimal_digits, output_format, row_labels, rel_uncertainty_threshold_labels, table_file_name, Print);
        movefile(table_file_name, output_folder);
    
        %--% Fuzzy points percentage table %--%
        mean_fuzzy_point_perc_list  = mean(fuzzy_points_percentage_matrix_b, 1);
        std_fuzzy_point_perc_list   = std(fuzzy_points_percentage_matrix_b, [], 1);
        max_fuzzy_point_perc_list   = max(fuzzy_points_percentage_matrix_b, [], 1);
        min_fuzzy_point_perc_list   = min(fuzzy_points_percentage_matrix_b, [], 1);
    
        fuzzy_point_perc_table_matrix   = [mean_fuzzy_point_perc_list; std_fuzzy_point_perc_list; max_fuzzy_point_perc_list; min_fuzzy_point_perc_list];
        row_labels                      = {'mean', 'std', 'max', 'min'};
    
        table_file_name = sprintf('Rel_Uncertainty_Fuzzy_Points_Percentage_%s.xls', bool_label);
        
        fprintf('Fuzzy points percentage for %s trees \n', bool_label);
        Table_Formatter(fuzzy_point_perc_table_matrix, num_decimal_digits, output_format, row_labels, rel_uncertainty_threshold_labels, table_file_name, Print);
        movefile(table_file_name, output_folder);
    
        %--% Combined plot %--%
        figure_number   = figure_number + 1;
        Fig             = figure(figure_number);
        figure_name     = sprintf('Percentages_Fuzzy_Cylinders_%s', bool_label);
    
        set(Fig, 'name', figure_name, 'NumberTitle', 'off');
    
        % Set the size and white background color
        set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85]);
        set(gcf, 'color', [1, 1, 1]);
        
        % Percentage of fuzzy cylinders
        subplot(1, 2, 1);
        hold on
        grid on
    
        plot(rel_uncertainty_threshold_steps, mean_fuzzy_cyl_perc_list, 'LineWidth', 2, 'color', 'b', 'DisplayName', 'mean');
    
        % sigma_lower_list = mean_fuzzy_cyl_perc_list - std_fuzzy_cyl_perc_list;
        % sigma_upper_list = mean_fuzzy_cyl_perc_list + std_fuzzy_cyl_perc_list;
        % patch('XData', [sigma_lower_list, fliplr(sigma_upper_list)]', 'YData', [rel_uncertainty_threshold_steps, fliplr(rel_uncertainty_threshold_steps)]', 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.5, 'DisplayName', sprintf('+/- 1%s', '\sigma'));
    
        xlim([0, max_threshold]);
        xlabel(sprintf('%s [-]', '\nu'));
        ylabel('Fuzzy fitted cylinders [%]');
    
        set(gca, 'FontSize', 15);
        set(gca, 'LineWidth', 2);
    
        hold off
    
        % Percentage of points for the fuzzy cylinders
        subplot(1, 2, 2);
        hold on
        grid on
    
        plot(rel_uncertainty_threshold_steps, mean_fuzzy_point_perc_list, 'LineWidth', 2, 'color', 'r', 'DisplayName', 'mean');
    
        % sigma_lower_list = mean_fuzzy_point_perc_list - std_fuzzy_point_perc_list;
        % sigma_upper_list = mean_fuzzy_point_perc_list + std_fuzzy_point_perc_list;
        % patch('XData', [sigma_lower_list, fliplr(sigma_upper_list)]', 'YData', [rel_uncertainty_threshold_steps, fliplr(rel_uncertainty_threshold_steps)]', 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.5, 'DisplayName', sprintf('+/- 1%s', '\sigma'));
    
        xlim([0, max_threshold]);
        xlabel(sprintf('%s [-]', '\nu'));
        ylabel('Points for fuzzy fitting [%]');
    
        set(gca, 'FontSize', 15);
        set(gca, 'LineWidth', 2);
    
        hold off
    
        % Saving the figure     
        export_fig(figure_number, sprintf('%s.fig', figure_name));
        export_fig(figure_number, sprintf('%s.png', figure_name));
        movefile([figure_name, '*'], output_folder); 
    end
    
%% Tree metric errors %%
    % Analysed tree metrics
    tree_metrics        = Results_Tables.Tree_Metrics.Rel_Error.Properties.VariableNames;
    number_tree_metrics = length(tree_metrics);

    for m = 1 : number_tree_metrics
        % This metric's data
        tree_metric = tree_metrics{m};

        if contains(tree_metric, 'Volume')
            % Volumes
            tree_metric_label = regexprep(tree_metric, '([A-Z])', ' ${lower($1)}');         % Inserts a space before each capital letter and makes them lower case

            if strcmp(tree_metric_label(1), ' ')                                            % This may result in a space being placed at the start, which is removed
                tree_metric_label(1) = [];
            end
    
            tree_metric_label = strrep(tree_metric_label, 'trunk', 'stem');                 % Trunk is replaced by stem
        else
            % DBH or height
            if strcmp(tree_metric, 'DBHqsm')
                tree_metric_label = 'DBH';
            elseif strcmp(tree_metric, 'UnmodDBHqsm')
                tree_metric_label = 'unmod. DBH';
            elseif strcmp(tree_metric, 'TreeHeight')
                tree_metric_label = 'tree height';
            end
        end

        Discr_Tree_Metrics_fun              = @(Results_Tables) str2double(Results_Tables.Tree_Metrics.Rel_Error('Discrete', tree_metric).Variables);
        discr_tree_metrics_rel_error_list   = cellfun(Discr_Tree_Metrics_fun, Results_Tables_cell(:, 1), 'UniformOutput', true);

        Fuzzy_Tree_Metrics_fun              = @(Results_Tables) str2double(Results_Tables.Tree_Metrics.Rel_Error('Fuzzy', tree_metric).Variables);
        fuzzy_tree_metrics_rel_error_matrix = cellfun(Fuzzy_Tree_Metrics_fun, Results_Tables_cell, 'UniformOutput', true);

        for b = 1 : number_bools
            % This boolean's data
            bool_list   = bool_cell{b};
            bool_label  = bool_label_cell{b};

            discr_tree_metrics_rel_error_list_b = discr_tree_metrics_rel_error_list(bool_list);
            discrete_mean_rel_error             = mean(discr_tree_metrics_rel_error_list_b);
    
            fuzzy_tree_metrics_rel_error_matrix_b   = fuzzy_tree_metrics_rel_error_matrix(bool_list, :);   
            fuzzy_mean_rel_error_list               = mean(fuzzy_tree_metrics_rel_error_matrix_b, 1);
            fuzzy_std_rel_error_list                = std(fuzzy_tree_metrics_rel_error_matrix_b, [], 1);
    
            %--% Figure %--%
            figure_number   = figure_number + 1;
            Fig             = figure(figure_number);
            figure_name     = sprintf('Relative_Error_%s_%s', tree_metric, bool_label);
      
            set(Fig, 'name', figure_name, 'NumberTitle', 'off');
    
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85]);
            set(gcf, 'color', [1, 1, 1]);
    
            hold on
            grid on
            
            % Discrete error
            plot(rel_uncertainty_threshold_steps, discrete_mean_rel_error * ones(1, number_steps), 'LineWidth', 2, 'color', 'b', 'DisplayName', 'TreeQSM mean');
    
            % Fuzzy error
            plot(rel_uncertainty_threshold_steps, fuzzy_mean_rel_error_list, 'LineWidth', 2, 'color', 'r', 'DisplayName', 'MixedQSM mean');
        
            % sigma_lower_list = fuzzy_mean_rel_error_list - fuzzy_std_rel_error_list;
            % sigma_upper_list = fuzzy_mean_rel_error_list + fuzzy_std_rel_error_list;
            % patch('XData', [rel_uncertainty_threshold_steps, fliplr(rel_uncertainty_threshold_steps)]', 'YData', [sigma_lower_list, fliplr(sigma_upper_list)]', 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.5, 'DisplayName', sprintf('+/- 1%s', '\sigma'));
    
            % Axes
            xlim([0, max_threshold]);
		    xlabel(sprintf('%s [-]', '\nu'));
            ylabel(sprintf('%s rel. error [%%]', tree_metric_label));
    
            % Legend
            legend('show', 'location', 'eastoutside');
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
    
            hold off
    
            % Saving the figure     
            export_fig(figure_number, sprintf('%s.fig', figure_name));
            export_fig(figure_number, sprintf('%s.png', figure_name));
            movefile([figure_name, '*'], output_folder); 
        end
    end

%% Branch order metrics %%
    % Metrics
    branch_order_metrics    = fieldnames(Results_Tables.Branch_Order_Metrics);
    number_metrics          = length(branch_order_metrics);

    % Colour map
    colour_map = cbrewer('seq', 'Reds', max(number_steps, 3));
    colour_map = max(0, colour_map);
    colour_map = min(1, colour_map);

    for m = 1 : number_metrics
        % This metric's data
        branch_order_metric = branch_order_metrics{m};

        Discr_BO_Metrics_fun                = @(Results_Tables) str2double(Results_Tables.Branch_Order_Metrics.(branch_order_metric).("Discr. rel. error"));
        discr_BO_metrics_rel_error_cell     = cellfun(Discr_BO_Metrics_fun, Results_Tables_cell(:, 1), 'UniformOutput', false);
        discr_BO_metrics_rel_error_matrix   = Padded_Vector_Matrix(discr_BO_metrics_rel_error_cell, NaN);

        Fuzzy_BO_Metrics_fun                = @(Results_Tables) str2double(Results_Tables.Branch_Order_Metrics.(branch_order_metric).("Fuzzy rel. error"));
        fuzzy_BO_metrics_rel_error_cell     = cellfun(Fuzzy_BO_Metrics_fun, Results_Tables_cell, 'UniformOutput', false);

        for b = 1 : number_bools
            % This boolean's data
            bool_list   = bool_cell{b};
            bool_label  = bool_label_cell{b};

            discr_BO_metrics_rel_error_matrix_b = discr_BO_metrics_rel_error_matrix(:, bool_list);
            discrete_mean_rel_error_list        = mean(discr_BO_metrics_rel_error_matrix_b, 2, 'omitnan');
            max_branch_order                    = length(discrete_mean_rel_error_list) - 1;
    
            fuzzy_mean_rel_error_matrix = zeros(max_branch_order + 1, number_steps);
            
            for i = 1 : number_steps
                fuzzy_BO_metrics_rel_error_cell_i   = fuzzy_BO_metrics_rel_error_cell(:, i);
                fuzzy_BO_metrics_rel_error_matrix_i = Padded_Vector_Matrix(fuzzy_BO_metrics_rel_error_cell_i, NaN);
                
                fuzzy_BO_metrics_rel_error_matrix_i(isinf(fuzzy_BO_metrics_rel_error_matrix_i)) = NaN;

                fuzzy_mean_rel_error_matrix(:, i)   = mean(fuzzy_BO_metrics_rel_error_matrix_i, 2, 'omitnan');
            end

            %--% Figure %--%
            figure_number   = figure_number + 1;
            Fig             = figure(figure_number);
            figure_name     = sprintf('Relative_Error_%s_%s', branch_order_metric, bool_label);
      
            set(Fig, 'name', figure_name, 'NumberTitle', 'off');
    
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85]);
            set(gcf, 'color', [1, 1, 1]);
    
            hold on
            grid on
            
            % Discrete error
            plot(0 : max_branch_order, discrete_mean_rel_error_list, 'LineWidth', 2, 'color', 'b', 'DisplayName', 'TreeQSM');
    
            % Fuzzy error
            for i = 1 : number_steps
                uncertainty_threshold       = rel_uncertainty_threshold_steps(i);
                colour                      = colour_map(i, :);
                fuzzy_mean_rel_error_list   = fuzzy_mean_rel_error_matrix(:, i);
    
                plot(0 : max_branch_order, fuzzy_mean_rel_error_list, 'LineWidth', 2, 'color', colour, 'DisplayName', sprintf('MixedQSM, %s = %.2g', '\nu', uncertainty_threshold));
            end
    
            % Axes
            xlim([0, max_branch_order]);
		    xlabel('branch order [-]');
            ylabel(sprintf('%s rel. error [%%]', branch_order_metric));
    
            % Legend
            legend('show', 'location', 'eastoutside');
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
    
            hold off
    
            % Saving the figure     
            export_fig(figure_number, sprintf('%s.fig', figure_name));
            export_fig(figure_number, sprintf('%s.png', figure_name));
            movefile([figure_name, '*'], output_folder); 
        end
    end
