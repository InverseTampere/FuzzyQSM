% This script tests fuzzy fitting the stem with fewer distributions to see if the extra computational time is worth it
% Note that it presumes that the TreeQSM QSMs have already been fitted

clear variables
close all
beep off
clc

%% Inputs %%
    %--% Sampling %--%
    max_number_distributions    =   1e3;                        % [-] Maximum number of distributions that's tested - excluding unlimited
    min_number_distributions    =   020;                        % [-] Minimum number of distributions that's tested
    number_steps                =   009;                        % [-] Number of sampling steps that are tested - including unlimited

    %--% Data %--%
    smoothed_models_folder      =   'Data/Smoothed_Models';                             % Folder in which the smoothed models are located
    TreeQSM_folder              =   'Data/TreeQSMs/030_beamdivergence';                 % Folder in which the TreeQSM generated QSMs are located
    point_clouds_folder         =   'Data/Point_Clouds/030_beamdivergence';             % Folder in which the point clouds are located
    output_folder               =   'Data/Results/030_beamdivergence/Stem_Fitting';     % Folder in which the results will be saved

    %--% Fuzzy fitting parameters %--%
    Fuzzy_Vector                =   true;                   % [true, false] Whether or not the vector is fitted or just taken from the discrete fit
    bounds_margin               =   0.50;                   % [-] Factor by which the geometry parameters can deviate from the initial estimates

    %--% Execution %--% 
    Parallel_Loop               =   true;                   % [true, false] Determines whether or not the Monte Carlo loop is ran in parallel
    max_cores                   =   10;                     % [-] If empty, the computer's number of cores is used. Otherwise this integer
    idle_timeout                =   030;                    % [min] Time until the parallel pool shuts down if idle        

    %--% Outputs %--%
    QSM_Print                   =   false;                  % [true, false] Shows printed statements regarding intermediate results
    QSM_Plot                    =   false;                  % [true, false] Shows plots of intermediate results
    Print                       =   true;                  % [true, false] Prints final results
    Plot                        =   true;                  % [true, false] Plots final results

    num_decimal_digits          =   4;                      % [-]
    output_format               =  'Float';                 % [Float, Integer, Exponential]
    
%% Input structures %%
    %--% Creation of the structures %--%
    Parallel_Pool           = struct('Parallel_Loop', Parallel_Loop, 'max_cores', max_cores, 'idle_timeout', idle_timeout);
    Fitting_Parameters      = struct('Fuzzy_Vector', Fuzzy_Vector, 'bounds_margin', bounds_margin);
    Output_Decisions        = struct('Print', QSM_Print, 'Plot', QSM_Plot, 'Compute_Radius_Errors', false);

%% Analysed trees %%
    % Tree IDs
    point_cloud_files   = dir(point_clouds_folder);
    point_cloud_files   = {point_cloud_files(:).name};

    dot_folders         = startsWith(point_cloud_files, '.');       % Remove the . and .. folders
    point_cloud_files   = point_cloud_files(~dot_folders);

    Tree_ID_fun     = @(point_cloud_file) strrep(point_cloud_file, '_Point_Cloud.mat', '');
    tree_ID_cell    = cellfun(Tree_ID_fun, point_cloud_files, 'UniformOutput', false);

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

%% Stem fitting %%
    % Maximum number of distributions
    max_number_distributions_list   = round(linspace(max_number_distributions, min_number_distributions, number_steps - 1));
    max_number_distributions_list   = [Inf, max_number_distributions_list];

    % Cycling through the trees
    Tree_Metrics_Rel_Errors_cell    = cell(number_trees, number_steps);
    comp_time_matrix                = zeros(number_trees, number_steps);
    percent_distr_matrix            = zeros(number_trees, number_steps);

    for t = 1 : number_trees
        % This tree's ID
        tree_ID = tree_ID_cell{t};

        % Its point cloud
        point_cloud_file_name   = point_cloud_files{t};
        Point_Cloud_File        = load(sprintf('%s/%s', point_clouds_folder, point_cloud_file_name));

        Point_Cloud_Data                        = Point_Cloud_File.Point_Cloud_Data_n;
        Scanning_Parameters                     = Point_Cloud_File.Scanning_Parameters;
        Scanner_Parameters                      = Point_Cloud_File.Scanner_Parameters;
        Scanner_Parameters.sigma_range_device   = Scanner_Parameters.sigma_range_0;

        % The TreeQSM fitted QSM
        TreeQSM_file_name   = sprintf('%s/%s_Discrete_QSM.mat', TreeQSM_folder, tree_ID);
        Tree_QSM_File       = load(TreeQSM_file_name);
        Discrete_QSM        = Tree_QSM_File.QSM_Discrete;
        TreeQSM_Inputs      = Tree_QSM_File.TreeQSM_Inputs;

        % The stem cylinders
        Discrete_Cylinders      = Discrete_QSM.cylinder;
        cyl_branch_order_list   = Discrete_QSM.cylinder.BranchOrder;
        stem_bool               = cyl_branch_order_list == 0;
        number_stem_cylinders   = sum(stem_bool);

        Stem_Cylinders          = Structure_Boolean(Discrete_Cylinders, stem_bool);
        Discrete_QSM.cylinder   = Stem_Cylinders;

        % Smoothed model metrics
        metrics_file_name           = sprintf('%s/%s/%s_Metrics.mat', smoothed_models_folder, tree_ID, tree_ID);
        Smoothed_Model_Metrics_File = load(metrics_file_name);
        Smoothed_Model_Metrics      = Smoothed_Model_Metrics_File.Smoothed_Model_Metrics;
        Original_QSM                = Smoothed_Model_Metrics_File.Original_QSM;
        Object_Circle_Geometry      = Smoothed_Model_Metrics_File.Object_Circle_Geometry;

        % The stem point clouds
        stem_point_indices_cell = Discrete_QSM.pmdistance.CylInd(stem_bool);
        number_points_list      = cellfun(@length, stem_point_indices_cell);
        number_points_total     = sum(number_points_list);

        % Fuzzy fitting for each number of distributions
        for i = 1 : number_steps
            % Reducing the number of points associated with each cylinder
            max_number_distributions    = max_number_distributions_list(i);
            stem_point_indices_cell_i   = stem_point_indices_cell;
    
            for c = 1 : number_stem_cylinders
                % Points of this cylinder
                point_indices = stem_point_indices_cell{c};
                number_points = number_points_list(c);
    
                if number_points > max_number_distributions
                    random_samples  = randperm(number_points, max_number_distributions);
                    point_indices   = point_indices(random_samples);
                    
                    stem_point_indices_cell_i{c} = point_indices;
                end
            end

            Discrete_QSM.pmdistance.CylInd = stem_point_indices_cell_i;

            % Percentage of remaining distributions
            number_points_list_i        = cellfun(@length, stem_point_indices_cell_i);
            number_points_total_i       = sum(number_points_list_i);
            percent_distr               = number_points_total_i / number_points_total * 100;
            percent_distr_matrix(t, i)  = percent_distr;
    
            % Fitting the fuzzy stem QSM
            tic;
    
            [Fuzzy_QSM_i, Cyl_Point_Cloud_Distributions_cell_i] = Fuzzy_QSM_Fitting(Discrete_QSM, Point_Cloud_Data, Fitting_Parameters, Scanner_Parameters, Scanning_Parameters, TreeQSM_Inputs, Parallel_Pool);
    
            comp_time               = toc;
            comp_time_matrix(t, i)  = comp_time;
    
            % Results
            Data_Parameters                     = struct();
            Results_Tables_i                    = FuzzyQSM_Result_Evaluation(tree_ID, Original_QSM, Discrete_QSM, Fuzzy_QSM_i, Smoothed_Model_Metrics, Object_Circle_Geometry, Cyl_Point_Cloud_Distributions_cell_i, Data_Parameters, Fitting_Parameters, Scanning_Parameters, Scanner_Parameters, Parallel_Pool, Output_Decisions);
            Tree_Metrics_Rel_Errors_cell{t, i}  = Results_Tables_i.Tree_Metrics.Rel_Error;
        end
    end

    % Save the results
    results_file_name = 'Tree_Metrics_Rel_Errors.mat';
    save(results_file_name, 'Tree_Metrics_Rel_Errors_cell', 'max_number_distributions_list', 'tree_ID_cell', 'comp_time_matrix', 'percent_distr_matrix');
    movefile(results_file_name, output_folder);

%% Results %%    
    %--% Computational time %--%
    % Relative computational time
    rel_comp_time_matrix    = comp_time_matrix ./ comp_time_matrix(:, 1) * 100;
    mean_rel_comp_time_list = mean(rel_comp_time_matrix, 1);

    % Plot
    figure_number = 1;
    Fig = figure(figure_number);
    
    figure_name = 'Num_Distr_Comp_Time';
    set(Fig, 'name', figure_name, 'NumberTitle', 'off');

    % Set the size and white background color
    set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
    set(gcf, 'color', [1, 1, 1])    
    
    hold on
    grid on

    plot(1:2, mean_rel_comp_time_list(1:2), 'LineWidth', 2, 'Color', 'b', 'LineStyle', ':');
    plot(2:number_steps, mean_rel_comp_time_list(2:number_steps), 'LineWidth', 2, 'Color', 'b');

    % Axes
    if number_steps ~= 1
        xlim([1, number_steps]);
        xticks(1 : number_steps);
        xlabel('max. nr. distributions [-]');
        xticklabels(['\infty', string(max_number_distributions_list(2 : number_steps))]);
    else
        xlabel('max nr. distributions');
    end
    
    ylabel('rel. computational time [%]');

    set(gca, 'FontSize', 15);
    set(gca, 'LineWidth', 2);
    
    hold off
    
    % Saving the figure    
    export_fig(figure_number, [figure_name, '.png']);
    movefile(sprintf('%s.png', figure_name), output_folder);

    %--% Percentage of remaining points/distributions %--%
    mean_percent_distr_list = mean(percent_distr_matrix, 1);

    % Plot
    figure_number = 2;
    Fig = figure(figure_number);
    
    figure_name = 'Num_Distr_Percent_Distr';
    set(Fig, 'name', figure_name, 'NumberTitle', 'off');

    % Set the size and white background color
    set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
    set(gcf, 'color', [1, 1, 1])    
    
    hold on
    grid on

    plot(1:2, mean_percent_distr_list(1:2), 'LineWidth', 2, 'Color', 'b', 'LineStyle', ':');
    plot(2:number_steps, mean_percent_distr_list(2:number_steps), 'LineWidth', 2, 'Color', 'b');

    % Axes
    if number_steps ~= 1
        xlim([1, number_steps]);
        xticks(1 : number_steps);
        xlabel('max. nr. distributions [-]');
        xticklabels(['\infty', string(max_number_distributions_list(2 : number_steps))]);
    else
        xlabel('max nr. distributions');
    end
    
    ylabel('remaining points [%]');

    set(gca, 'FontSize', 15);
    set(gca, 'LineWidth', 2);
    
    hold off
    
    % Saving the figure    
    export_fig(figure_number, [figure_name, '.png']);
    movefile(sprintf('%s.png', figure_name), output_folder);

    %--% The relative errors for the relevant tree metrics %--%
    tree_metrics    = {'TrunkVolume', 'DBHqsm'};           % The others don't make sense to analyse as no branches were fitted
    metric_labels   = {'Stem volume', 'DBH'};
    number_metrics  = length(tree_metrics);

    for m = 1 : number_metrics
        % This metric's data
        tree_metric     = tree_metrics{m};
        metric_label    = metric_labels{m};

        for b = 1 : number_bools
            % This boolean's data
            bool_list   = bool_cell{b};
            bool_label  = bool_label_cell{b};
            
            % Discrete relative errors for each tree
            Discr_Rel_Error_fun     = @(Tree_Metrics_Rel_Error) str2double(Tree_Metrics_Rel_Error('Discrete', tree_metric).Variables);
            discr_rel_error_list    = cellfun(Discr_Rel_Error_fun, Tree_Metrics_Rel_Errors_cell(bool_list, 1));                     
    
            [discr_mean, discr_max, discr_min, discr_std] = deal(mean(discr_rel_error_list), max(discr_rel_error_list), min(discr_rel_error_list), std(discr_rel_error_list));
    
            % Fuzzy relative errors for each tree and max. number of distributions
            Fuzzy_Rel_Error_fun     = @(Tree_Metrics_Rel_Error) str2double(Tree_Metrics_Rel_Error('Fuzzy', tree_metric).Variables);
            fuzzy_rel_error_matrix  = cellfun(Fuzzy_Rel_Error_fun, Tree_Metrics_Rel_Errors_cell(bool_list, :));
    
            [fuzzy_mean_list, fuzzy_max_list, fuzzy_min_list, fuzzy_std_list] = deal(mean(fuzzy_rel_error_matrix, 1), max(fuzzy_rel_error_matrix, [], 1), min(fuzzy_rel_error_matrix, [], 1), std(fuzzy_rel_error_matrix, [], 1));
    
            % Table
            Label_fun                       = @(n) sprintf('Fuzzy %i', n);
            max_number_distributions_labels = arrayfun(Label_fun, max_number_distributions_list, 'UniformOutput', false);
            row_labels                      = ['Discrete', max_number_distributions_labels];
            
            column_labels   = {'mean', 'max', 'min', 'std'};
            table_matrix    = [discr_mean, discr_max, discr_min, discr_std; 
                               fuzzy_mean_list', fuzzy_max_list', fuzzy_min_list', fuzzy_std_list'];
    
            if Fuzzy_Vector == true
                table_file_name = sprintf('%s_%s_Num_Distr_Rel_Errors_FuzVec.xls', tree_metric, bool_label);
            else
                table_file_name = sprintf('%s_%s_Num_Distr_Rel_Errors.xls', tree_metric, bool_label);
            end

            fprintf('Rel. error of %s for %s trees \n', tree_metric, bool_label);
            Table_Formatter(table_matrix, num_decimal_digits, output_format, row_labels, column_labels, table_file_name, Print);
            movefile(table_file_name, output_folder);
    
            % Plot
            figure_number = figure_number + 1;
            Fig = figure(figure_number);
            
            if Fuzzy_Vector == true
                figure_name = sprintf('%s_%s_Num_Distr_Rel_Errors_FuzVec.png', tree_metric, bool_label);
            else
                figure_name = sprintf('%s_%s_Num_Distr_Rel_Errors.png', tree_metric, bool_label);
            end

            set(Fig, 'name', figure_name, 'NumberTitle', 'off');
        
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            hold on
            grid on
            
            % Zero line for clarity
            pl_zero = plot([1, number_steps], [0, 0], 'LineWidth', 2, 'color', 'k', 'DisplayName', 'zero-line');
            pl_zero.HandleVisibility = 'Off';
    
            % Discrete mean rel. error 
            plot([0, 1], discr_mean * [1, 1], 'LineWidth', 2, 'Color', 'r', 'LineStyle', ':', 'HandleVisibility', 'Off');
            plot([1, number_steps], discr_mean * [1, 1], 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'TreeQSM mean');
    
            % Fuzzy mean and standard deviation
            plot(1 : 2, fuzzy_mean_list(1:2), 'LineWidth', 2, 'Color', 'b', 'LineStyle', ':', 'HandleVisibility', 'Off');
            plot(2 : number_steps, fuzzy_mean_list(2:end), 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'FuzzyQSM mean');
    
            % sigma_lower_list = fuzzy_mean_list - fuzzy_std_list;
            % sigma_upper_list = fuzzy_mean_list + fuzzy_std_list;
            % patch('XData', [1 : number_steps, number_steps : - 1 : 1]', 'YData', [sigma_lower_list, fliplr(sigma_upper_list)]', 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.5, 'DisplayName', sprintf('+/- 1%s', '\sigma'));
    
            % Axes
            if number_steps ~= 1
                xlim([1, number_steps]);
                xticks(1 : number_steps);
                xlabel('max. nr. distributions [-]');
                xticklabels(['\infty', string(max_number_distributions_list(2 : number_steps))]);
            else
                xlabel('max nr. distributions');
            end
        
            ylabel(sprintf('%s relative error [%%]', metric_label));
        
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
        
            hold off
        
            % Legend
            legend('show', 'location', 'eastoutside');
        
            % Saving the figure    
            export_fig(figure_number, [figure_name, '.png']);
            movefile(sprintf('%s.png', figure_name), output_folder);
        end
    end