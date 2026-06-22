% The results are evaluated for fuzzy fitting of the Wytham trees

clear variables
close all
beep off
clc

%% Inputs %%
    % Execution
    Parallel_Loop               =   false;                   % [true, false] Determines whether or not certain metrics are determined in parallel
    max_cores                   =   10;                     % [-] If empty, the computer's number of cores is used. Otherwise this integer
    idle_timeout                =   030;                    % [min] Time until the parallel pool shuts down if idle        

    % File locations
    discrete_QSMs_folder        =   'Data/TreeQSMs/030_beamdivergence';        % String of the folder in which the discrete QSMs are located
    fuzzy_QSMs_folder           =   'Data/FuzzyQSMs/030_beamdivergence';       % String of the folder in which the fuzzy QSMs are located
    smoothed_models_folder      =   'Data/Smoothed_Models';                     % Folder in which the smoothed models are located
    results_folder              =   'Data/Results/030_beamdivergence';         % Folder in which the results will be saved
        
    % Outputs
    Individual_Tree_Outputs     =   true;                  % [true, false] Shows results for each single tree
    Compute_Radius_Errors       =   false;                  % [true, false] Whether or not the radius errors are computed, which are slow
    num_decimal_digits          =   4;                      % [-]
    output_format               =   'Float';                % [Float, Integer, Exponential]
    Print                       =   true;                  % [true, false] Whether or not the tables are displayed

%% Analysed trees %%
    % Tree ID folders which contain the data
    smoothed_model_folders   = dir(smoothed_models_folder);
    smoothed_model_folders   = {smoothed_model_folders(:).name};

    dot_folders     = contains(smoothed_model_folders, '.');       % Remove the . and .. folders
    tree_ID_cell    = smoothed_model_folders(~dot_folders);
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

    % And for all categories combined
    bool_cell       = [{true(1, number_trees)}, species_bool_cell];
    bool_label_cell = ['All', unique_tree_species];
    number_bools    = length(bool_cell);

%% Individual tree results %%
    % Results for each tree
    Results_Tables_cell = cell(number_trees, 1);

    for t = 1 : number_trees
        %--% This tree's data %--%
        % Its ID (class and number)
        tree_ID         = tree_ID_cell{t};

        % Check if the results were already created
        tree_folder     = dir(sprintf('%s_*', tree_ID));
        number_folders  = length(tree_folder);

        if number_folders == 1
            % Load them if they are present
            results_file_name       = sprintf('%s/FuzzyQSM_Outputs.mat', tree_folder.name);
            Results_File            = load(results_file_name);
            Results_Tables          = Results_File.Results_Tables;
            Results_Tables_cell{t}  = Results_Tables;

            continue
        end

        % Discrete QSM
        discrete_QSM_file_name  = sprintf('%s/%s_Discrete_QSM.mat', discrete_QSMs_folder, tree_ID);
        Discrete_QSM_File       = load(discrete_QSM_file_name);
        Discrete_QSM            = Discrete_QSM_File.QSM_Discrete;

        % Fuzzy QSM
        fuzzy_QSM_file_name     = sprintf('%s/%s_Fuzzy_QSM.mat', fuzzy_QSMs_folder, tree_ID);
        Fuzzy_QSM_File          = load(fuzzy_QSM_file_name);
        Fuzzy_QSM               = Fuzzy_QSM_File.QSM_Fuzzy;

        % Smoothed model metrics
        metrics_file_name                   = sprintf('%s/%s/%s_Metrics.mat', smoothed_models_folder, tree_ID, tree_ID);
        Smoothed_Model_Metrics_File         = load(metrics_file_name);
        Smoothed_Model_Metrics              = Smoothed_Model_Metrics_File.Smoothed_Model_Metrics;
        Smoothed_Model_Circle_Geometry      = Smoothed_Model_Metrics_File.Smoothed_Model_Circle_Geometry;
        Original_QSM                        = Smoothed_Model_Metrics_File.Original_QSM;

        % Distributions
        distributions_file_name             = sprintf('%s/%s_Point_Cloud_Distributions.mat', fuzzy_QSMs_folder, tree_ID);
        Distributions_File                  = load(distributions_file_name);
        Cyl_Point_Cloud_Distributions_cell  = Distributions_File.Cyl_Point_Cloud_Distributions_cell;

        %--% Results for this tree %--%
        Output_Decisions        = struct('Print', Individual_Tree_Outputs, 'Plot', Individual_Tree_Outputs, 'Compute_Radius_Errors', Compute_Radius_Errors);
        Parallel_Pool           = struct('Parallel_Loop', Parallel_Loop, 'max_cores', max_cores, 'idle_timeout', idle_timeout);
        Results_Tables          = FuzzyQSM_Result_Evaluation(tree_ID, Original_QSM, Discrete_QSM, Fuzzy_QSM, Smoothed_Model_Metrics, Smoothed_Model_Circle_Geometry, Cyl_Point_Cloud_Distributions_cell, [], [], [], [], Parallel_Pool, Output_Decisions);
        Results_Tables_cell{t}  = Results_Tables;

        if Individual_Tree_Outputs == true
            close all
        end
    end

%% Combined tree metrics %%
    % Assessed tree metrics
    tree_metrics        = Results_Tables.Tree_Metrics.Rel_Error.Properties.VariableNames;
    number_tree_metrics = length(tree_metrics);

    %--% Relative error tables %--%
    % One table per boolean
    for b = 1 : number_bools
        % The boolean
        bool_list   = bool_cell{b};
        bool_label  = bool_label_cell{b};

        % Relative errors of the discrete and fuzzy QSMs
        Discr_Tree_Metrics_Rel_Error_fun    = @(Results_Tables) str2double(Results_Tables.Tree_Metrics.Rel_Error('Discrete', :).Variables);
        discr_tree_metrics_rel_error_cell   = cellfun(Discr_Tree_Metrics_Rel_Error_fun, Results_Tables_cell(bool_list), 'UniformOutput', false);
        discr_tree_metrics_rel_error_matrix = vertcat(discr_tree_metrics_rel_error_cell{:});
        
        Fuzzy_Tree_Metrics_Rel_Error_fun    = @(Results_Tables) str2double(Results_Tables.Tree_Metrics.Rel_Error('Fuzzy', :).Variables);
        fuzzy_tree_metrics_rel_error_cell   = cellfun(Fuzzy_Tree_Metrics_Rel_Error_fun, Results_Tables_cell(bool_list), 'UniformOutput', false);
        fuzzy_tree_metrics_rel_error_matrix = vertcat(fuzzy_tree_metrics_rel_error_cell{:});
    
        % The mean (absolute), maximum and minimum errors for each tree metric are tabulated
        table_row_labels    = {'Discr. mean error [%]', 'Discr. mean abs. error [%]', 'Discr. max. error [%]', 'Discr. min. error [%]', ...
                               'Fuzzy mean error [%]', 'Fuzzy mean abs. error [%]', 'Fuzzy max. error [%]', 'Fuzzy min. error [%]'};
        tree_metrics_data   = [mean(discr_tree_metrics_rel_error_matrix, 1); mean(abs(discr_tree_metrics_rel_error_matrix), 1); max(discr_tree_metrics_rel_error_matrix, [], 1); min(discr_tree_metrics_rel_error_matrix, [], 1); ...
                               mean(fuzzy_tree_metrics_rel_error_matrix, 1); mean(abs(fuzzy_tree_metrics_rel_error_matrix), 1); max(fuzzy_tree_metrics_rel_error_matrix, [], 1); min(fuzzy_tree_metrics_rel_error_matrix, [], 1)];
                
        table_file_name     = sprintf('Tree_Metrics_%s_Relative_Errors.xls', bool_label);
        
        fprintf('Tree metrics for %s trees: \n', bool_label);
        Table_Formatter(tree_metrics_data, num_decimal_digits, output_format, table_row_labels, tree_metrics, table_file_name, Print);
    
        movefile(table_file_name, results_folder);
    end

    % One table per metric
    for m = 1 : number_tree_metrics
        tree_metric = tree_metrics{m};

        discr_mean_rel_error_list       = zeros(1, number_bools);
        discr_mean_abs_rel_error_list   = zeros(1, number_bools);
        fuzzy_mean_rel_error_list       = zeros(1, number_bools);
        fuzzy_mean_abs_rel_error_list   = zeros(1, number_bools);

        for b = 1 : number_bools
            % The boolean
            bool_list   = bool_cell{b};

            % Relative errors of the discrete and fuzzy QSMs
            Discr_Tree_Metrics_Rel_Error_fun    = @(Results_Tables) str2double(Results_Tables.Tree_Metrics.Rel_Error('Discrete', tree_metric).Variables);
            discr_tree_metrics_rel_error_list   = cellfun(Discr_Tree_Metrics_Rel_Error_fun, Results_Tables_cell(bool_list), 'UniformOutput', true);
            discr_mean_rel_error_list(b)        = mean(discr_tree_metrics_rel_error_list);
            discr_mean_abs_rel_error_list(b)    = mean(abs(discr_tree_metrics_rel_error_list));

            Fuzzy_Tree_Metrics_Rel_Error_fun    = @(Results_Tables) str2double(Results_Tables.Tree_Metrics.Rel_Error('Fuzzy', tree_metric).Variables);
            fuzzy_tree_metrics_rel_error_list   = cellfun(Fuzzy_Tree_Metrics_Rel_Error_fun, Results_Tables_cell(bool_list), 'UniformOutput', true);
            fuzzy_mean_rel_error_list(b)        = mean(fuzzy_tree_metrics_rel_error_list);
            fuzzy_mean_abs_rel_error_list(b)    = mean(abs(fuzzy_tree_metrics_rel_error_list));
        end

        % The mean (absolute), maximum and minimum errors for each tree metric are tabulated
        table_row_labels    = {'Discrete mean rel. error [%]', 'Discrete mean abs. rel. error [%]', 'Fuzzy mean rel. error [%]', 'Fuzzy abs. mean rel. error [%]'};
        tree_metric_data    = [discr_mean_rel_error_list; discr_mean_abs_rel_error_list; fuzzy_mean_rel_error_list; fuzzy_mean_abs_rel_error_list];
                
        table_file_name     = sprintf('Tree_Metric_%s_Relative_Errors.xls', tree_metric);
        
        fprintf('%s: \n', tree_metric);
        Table_Formatter(tree_metric_data, num_decimal_digits, output_format, table_row_labels, bool_label_cell, table_file_name, Print);
    
        movefile(table_file_name, results_folder);
    end

    %--% Value plots %--%
    % Colour map for each species
    species_cmap = cbrewer('qual', 'Set1', max(3, number_species));

    % Tree metrics data
    QSM_labels          = {'Discrete', 'Fuzzy', 'Original'};
    number_QSMs         = length(QSM_labels);
    tree_metrics_cell   = cell(1, number_QSMs);

    for q = 1 : number_QSMs
        % This QSM's tree metrics data
        QSM_label               = QSM_labels{q};
        Tree_Metrics_fun        = @(Results_Tables) cellfun(@str2double, Results_Tables.Tree_Metrics.Values{QSM_label, :});
        tree_metrics_cell_q     = cellfun(Tree_Metrics_fun, Results_Tables_cell, 'UniformOutput', false);
        tree_metrics_matrix_q   = vertcat(tree_metrics_cell_q{:});
        tree_metrics_cell{q}    = tree_metrics_matrix_q;
    end

    [discrete_tree_metrics_matrix, fuzzy_tree_metrics_matrix, original_tree_metrics_matrix] = deal(tree_metrics_cell{:});

    for t = 1 : number_tree_metrics
        % This metric's data
        tree_metric                 = tree_metrics{t};
        discrete_tree_metrics_list  = discrete_tree_metrics_matrix(:, t);
        fuzzy_tree_metrics_list     = fuzzy_tree_metrics_matrix(:, t);
        original_tree_metrics_list  = original_tree_metrics_matrix(:, t);
        
        % Data bounds with margin
        min_value   = min([discrete_tree_metrics_list; fuzzy_tree_metrics_list; original_tree_metrics_list]);
        max_value   = max([discrete_tree_metrics_list; fuzzy_tree_metrics_list; original_tree_metrics_list]);
        data_ampl   = max_value - min_value;
        data_LB     = min_value - 0.1*data_ampl;
        data_UB     = max_value + 0.1*data_ampl;

        % Neat label and unit
        if contains(tree_metric, 'Volume')
            % Volumes
            tree_metric_unit = 'm^3';

            tree_metric_label = regexprep(tree_metric, '([A-Z])', ' ${lower($1)}');         % Inserts a space before each capital letter and makes them lower case

            if strcmp(tree_metric_label(1), ' ')                                            % This may result in a space being placed at the start, which is removed
                tree_metric_label(1) = [];
            end
    
            tree_metric_label = strrep(tree_metric_label, 'trunk', 'stem');                 % Trunk is replaced by stem
        else
            % DBH, height or length
            tree_metric_unit = 'm';

            if strcmp(tree_metric, 'DBHqsm')
                tree_metric_label = 'DBH';
            elseif strcmp(tree_metric, 'UnmodDBHqsm')
                tree_metric_label = 'unmod. DBH';
            elseif strcmp(tree_metric, 'TreeHeight')
                tree_metric_label = 'tree height';
            elseif strcmp(tree_metric, 'BranchLength')
                tree_metric_label = 'branch length';
            end
        end

        %--% The figure %--%
        figure_number   = t;
        Fig             = figure(figure_number);
        figure_name     = tree_metric_label;
    
        set(Fig, 'name', figure_name, 'NumberTitle', 'off');

        % Set the size and white background color
        set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
        set(gcf, 'color', [1, 1, 1])    

        hold on
        grid on

        % Fuzzy and discrete legend entries
        scatter(NaN, NaN, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', 'TreeQSM');
        scatter(NaN, NaN, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', 'FuzzyQSM');

        % Fuzzy and discrete values vs. the original values
        for s = 1 : number_species
            % This species' data
            species_colour  = species_cmap(s, :);
            species_bool    = species_bool_cell{s};
            species_label   = unique_tree_species{s};

            % Scatter plots
            scatter(original_tree_metrics_list(species_bool), discrete_tree_metrics_list(species_bool), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', species_colour, 'LineWidth', 2, 'HandleVisibility', 'Off');
            scatter(original_tree_metrics_list(species_bool), fuzzy_tree_metrics_list(species_bool), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', species_colour, 'LineWidth', 2, 'HandleVisibility', 'Off');

            % Marker for the legend
            scatter(NaN, NaN, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', species_colour, 'LineWidth', 2, 'DisplayName', species_label);
        end

        % Diagonal
        pl_diag = plot([data_LB, data_UB], [data_LB, data_UB], 'LineWidth', 2, 'LineStyle', ':');
        pl_diag.HandleVisibility = 'Off';

        % Axes
        xlim([data_LB, data_UB]);
        ylim([data_LB, data_UB]);

        xlabel(sprintf('Original %s [%s]', tree_metric_label, tree_metric_unit));
        ylabel(sprintf('Fitted %s [%s]', tree_metric_label, tree_metric_unit));

        % Legend
        legend('show', 'location', 'eastoutside');

        set(gca, 'FontSize', 15);
        set(gca, 'LineWidth', 2);

        hold off

        % Saving the figure        
        export_fig(figure_number, sprintf('%s.fig', figure_name));
        export_fig(figure_number, sprintf('%s.png', figure_name));

        movefile([figure_name, '*'], results_folder); 
    end

%% Combined branch order metrics %%
    % Assessed branch order metrics
    branch_order_metrics    = fieldnames(Results_Tables.Branch_Order_Metrics);
    number_metrics          = length(branch_order_metrics);

    % Maximum branch order
    branch_order_metric     = branch_order_metrics{1};
    Max_Branch_Order_fun    = @(Results_Tables) str2double(Results_Tables.Branch_Order_Metrics.(branch_order_metric).Row{end});
    max_branch_order_list   = cellfun(Max_Branch_Order_fun, Results_Tables_cell);
    max_branch_order        = max(max_branch_order_list);
    branch_order_labels     = string(0 : max_branch_order);

    for m = 1 : number_metrics
        % This metric's errors
        branch_order_metric = branch_order_metrics{m};

        discr_metric_rel_error_matrix = NaN(max_branch_order + 1, number_trees);
        fuzzy_metric_rel_error_matrix = NaN(max_branch_order + 1, number_trees);

        for t = 1 : number_trees
            % This tree's results
            Results_Tables = Results_Tables_cell{t};

            % Discrete and fuzzy errors
            discr_metric_rel_error_cell = Results_Tables.Branch_Order_Metrics.(branch_order_metric).("Discr. rel. error");
            discr_metric_rel_error_list = cellfun(@str2double, discr_metric_rel_error_cell);

            fuzzy_metric_rel_error_cell = Results_Tables.Branch_Order_Metrics.(branch_order_metric).("Fuzzy rel. error");
            fuzzy_metric_rel_error_list = cellfun(@str2double, fuzzy_metric_rel_error_cell);

            number_branch_orders                                        = length(discr_metric_rel_error_list);
            discr_metric_rel_error_matrix(1 : number_branch_orders, t)  = discr_metric_rel_error_list;
            fuzzy_metric_rel_error_matrix(1 : number_branch_orders, t)  = fuzzy_metric_rel_error_list;
        end

        % Infinites are replaced by NaN
        discr_metric_rel_error_matrix(isinf(discr_metric_rel_error_matrix)) = NaN;
        fuzzy_metric_rel_error_matrix(isinf(fuzzy_metric_rel_error_matrix)) = NaN;

        for b = 1 : number_bools
            % The boolean's data
            bool_list   = bool_cell{b};
            bool_label  = bool_label_cell{b};

            discr_metric_rel_error_matrix_b = discr_metric_rel_error_matrix(:, bool_list);
            fuzzy_metric_rel_error_matrix_b = fuzzy_metric_rel_error_matrix(:, bool_list);

            % The mean (absolute), maximum and minimum errors for each branch order metric metric are tabulated
            table_row_labels            = {'Discr. mean error [%]', 'Discr. mean abs. error [%]', 'Discr. max. error [%]', 'Discr. min. error [%]', ...
                                           'Fuzzy mean error [%]', 'Fuzzy mean abs. error [%]', 'Fuzzy max. error [%]', 'Fuzzy min. error [%]'};
            branch_order_metrics_data   = [mean(discr_metric_rel_error_matrix_b, 2, 'omitnan'), mean(abs(discr_metric_rel_error_matrix_b), 2, 'omitnan'), max(discr_metric_rel_error_matrix_b, [], 2), min(discr_metric_rel_error_matrix_b, [], 2), ...
                                           mean(fuzzy_metric_rel_error_matrix_b, 2, 'omitnan'), mean(abs(fuzzy_metric_rel_error_matrix_b), 2, 'omitnan'), max(fuzzy_metric_rel_error_matrix_b, [], 2), min(fuzzy_metric_rel_error_matrix_b, [], 2)]';
                                
            table_file_name             = sprintf('BO_%s_%s_Relative_Errors.xls', branch_order_metric, bool_label);

            fprintf('%s for %s trees: \n', branch_order_metric, bool_label);
            Table_Formatter(branch_order_metrics_data, num_decimal_digits, output_format, table_row_labels, branch_order_labels, table_file_name, Print);
        
            movefile(table_file_name, results_folder);
        end
    end

%% Combined radius errors %%
    if Compute_Radius_Errors == true
        % Assessed radius error metrics
        radius_error_metrics    = fieldnames(Results_Tables.Radius);
        number_errors           = length(radius_error_metrics);
    
        for e = 1 : number_errors
            % Error metric
            error_metric = radius_error_metrics{e};
    
            % Errors for each tree
            discr_radius_error_matrix = NaN(max_branch_order + 2, number_trees);            % Note that the total is included
            fuzzy_radius_error_matrix = NaN(max_branch_order + 2, number_trees);
        
            for t = 1 : number_trees
                % This tree's results
                Results_Tables = Results_Tables_cell{t};
        
                % Discrete and fuzzy errors
                discr_radius_error_cell = Results_Tables.Radius.(error_metric).(sprintf('Discrete mean radius %s', error_metric));
                discr_radius_error_list = cellfun(@str2double, discr_radius_error_cell);
        
                fuzzy_radius_error_cell = Results_Tables.Radius.(error_metric).(sprintf('Fuzzy mean radius %s', error_metric));
                fuzzy_radius_error_list = cellfun(@str2double, fuzzy_radius_error_cell);
        
                number_branch_orders                                    = length(discr_radius_error_list);
                discr_radius_error_matrix(1 : number_branch_orders, t)  = discr_radius_error_list;
                fuzzy_radius_error_matrix(1 : number_branch_orders, t)  = fuzzy_radius_error_list;
            end
        
            % The mean (absolute), maximum and minimum errors for each tree metric are tabulated
            for b = 1 : number_bools
                % This boolean's data
                bool_list   = bool_cell{b};
                bool_label  = bool_label_cell{b};
    
                discr_radius_error_matrix_b = discr_radius_error_matrix(:, bool_list);
                fuzzy_radius_error_matrix_b = fuzzy_radius_error_matrix(:, bool_list);
    
                % The mean (absolute), maximum and minimum errors for each radius error metric are tabulated
                table_column_labels = ['Total', branch_order_labels];

                table_row_labels    = {'Discr. mean error [%]', 'Discr. mean abs. error [%]', 'Discr. max. error [%]', 'Discr. min. error [%]', ...
                                       'Fuzzy mean error [%]', 'Fuzzy mean abs. error [%]', 'Fuzzy max. error [%]', 'Fuzzy min. error [%]'};
                radius_error_data   = [mean(discr_radius_error_matrix_b, 2, 'omitnan'), mean(abs(discr_radius_error_matrix_b), 2, 'omitnan'), max(discr_radius_error_matrix_b, [], 2), min(discr_radius_error_matrix_b, [], 2), ...
                                       mean(fuzzy_radius_error_matrix_b, 2, 'omitnan'), mean(abs(fuzzy_radius_error_matrix_b), 2, 'omitnan'), max(fuzzy_radius_error_matrix_b, [], 2), min(fuzzy_radius_error_matrix_b, [], 2)]';
                                    
                table_file_name     = sprintf('Radius_%s_%s.xls', error_metric, bool_label);
    
                fprintf('%s radius error for %s trees: \n', error_metric, bool_label);
                Table_Formatter(radius_error_data, num_decimal_digits, output_format, table_row_labels, table_column_labels, table_file_name, Print);
            
                movefile(table_file_name, results_folder);
            end
        end
    end

%% Combined stem taper curve %%
    if Compute_Radius_Errors == true
        % Analysed stem taper error metrics
        stem_error_metrics  = Results_Tables.Stem_Taper.Properties.VariableNames;
        number_stem_metrics = length(stem_error_metrics);
    
        for b = 1 : number_bools
            % This boolean's data
            bool_list   = bool_cell{b};
            bool_label  = bool_label_cell{b};
    
            % The stem errors
            discrete_stem_error_list    = zeros(1, number_stem_metrics);
            fuzzy_stem_error_list       = zeros(1, number_stem_metrics);
        
            for s = 1 : number_stem_metrics
                % This stem's error metric
                stem_metric = stem_error_metrics{s};
        
                % Discrete value
                Discrete_Stem_Error_fun     = @(Results_Tables) Results_Tables.Stem_Taper('Discrete', stem_metric).Variables;
                discrete_stem_error_cell_s  = cellfun(Discrete_Stem_Error_fun, Results_Tables_cell(bool_list), 'UniformOutput', false);
                discrete_stem_error_list_s  = cellfun(@str2double, discrete_stem_error_cell_s);
                discrete_mean_stem_error_s  = mean(discrete_stem_error_list_s);
                discrete_stem_error_list(s) = discrete_mean_stem_error_s;
        
                % Fuzzy value
                Fuzzy_Stem_Error_fun        = @(Results_Tables) Results_Tables.Stem_Taper('Fuzzy', stem_metric).Variables;
                fuzzy_stem_error_cell_s     = cellfun(Fuzzy_Stem_Error_fun, Results_Tables_cell(bool_list), 'UniformOutput', false);
                fuzzy_stem_error_list_s     = cellfun(@str2double, fuzzy_stem_error_cell_s);
                fuzzy_mean_stem_error_s     = mean(fuzzy_stem_error_list_s);
                fuzzy_stem_error_list(s)    = fuzzy_mean_stem_error_s;
            end
            
            % Table
            table_row_labels    = {sprintf('TreeQSM %s', bool_label), sprintf('FuzzyQSM %s', bool_label)};
            stem_table_matrix   = [discrete_stem_error_list; fuzzy_stem_error_list] * 1e3;      % m to mm

            m_to_mm_fun             = @(metric) strrep(metric, '[m]', '[mm]');
            stem_error_metrics_mm   = cellfun(m_to_mm_fun, stem_error_metrics, 'UniformOutput', false);

            table_file_name     = sprintf('Stem_Taper_Errors_%s.xls', bool_label);
            Table_Formatter(stem_table_matrix, num_decimal_digits, output_format, table_row_labels, stem_error_metrics_mm, table_file_name, Print);
        
            movefile(table_file_name, results_folder);
        end
    end