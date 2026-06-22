% This script determines the correlation between a given metric of the original QSMs and absolute relative errors

clear variables
close all
beep off
clc

%% Inputs %%
    % Correlation metric
    correlation_metric          = 'BranchLength';                 % Name of the tree metric field in the original QSM against which is tested

    % Beam divergence folders
    beam_divergence_1_folder    = 'Data/Results/030_beamdivergence';
    beam_divergence_2_folder    = 'Data/Results/0175_beamdivergence';

    beam_divergence_colours     = {'b', 'r'};

    % Outputs
    output_folder               = 'Data/Results/Error_Correlation/BranchLength';
    num_decimal_digits          = 4;                        % [-]
    output_format               = 'Float';                  % [Float, Integer, Exponential]
    Print                       = true;                     % [true, false] Whether or not the tables are displayed

%% Tree metrics %%
    % Beam divergences
    beam_divergence_folders     = {beam_divergence_1_folder, beam_divergence_2_folder};
    number_beam_divergences     = length(beam_divergence_folders);

    beam_divergence_list = zeros(1, number_beam_divergences);

    for b = 1 : number_beam_divergences
        folder_name     = beam_divergence_folders{b};
        name_parts      = strsplit(folder_name, '/');

        Beam_Div_Part_fun   = @(part) contains(part, 'beamdivergence');
        beam_div_bool       = cellfun(Beam_Div_Part_fun, name_parts);
        beam_div_part       = name_parts{beam_div_bool};

        beam_divergence = strrep(beam_div_part, '_beamdivergence', '');
        beam_divergence = sprintf('%s.%s', beam_divergence(1), beam_divergence(2 : end));
        beam_divergence = str2double(beam_divergence);

        beam_divergence_list(b) = beam_divergence;
    end

    [beam_divergence_list, order]   = sort(beam_divergence_list);
    beam_divergence_folders         = beam_divergence_folders(order);

    % Tree metrics
    Tree_Metrics_cell = cell(1, number_beam_divergences);

    for b = 1 : number_beam_divergences
        beam_divergence_folder  = beam_divergence_folders{b};
        tree_folders            = dir(beam_divergence_folder);
        tree_folders            = {tree_folders(:).name};
    
        dot_folders             = contains(tree_folders, '.');       % Remove the . and .. folders and any files
        tree_folders            = tree_folders(~dot_folders);

        relative_uncertainty_bool   = contains(tree_folders, 'Relative_Uncertainty');
        tree_folders                = tree_folders(~relative_uncertainty_bool);

        stem_fitting_bool           = contains(tree_folders, 'Stem_Fitting');
        tree_folders                = tree_folders(~stem_fitting_bool);

        number_trees            = length(tree_folders);

        Tree_Metrics_cell_b = cell(1, number_trees);

        for t = 1 : number_trees
            % This tree's metrics
            tree_folder             = tree_folders{t};
            tree_results_file       = sprintf('%s/%s/FuzzyQSM_Outputs.mat', beam_divergence_folder, tree_folder);
            Tree_Results            = load(tree_results_file);
    
            Tree_Metrics            = Tree_Results.Results_Tables.Tree_Metrics;
            Tree_Metrics_cell_b{t}  = Tree_Metrics;
        end

        Tree_Metrics_cell{b} = Tree_Metrics_cell_b;
    end

    % Change the array structure
    Tree_Metrics_cell = vertcat(Tree_Metrics_cell{:});

%% Correlation %%
    % Original correlation metric values
    Original_Corr_Metric_fun    = @(Tree_Metrics) str2double(Tree_Metrics.Values('Original', correlation_metric).Variables);
    original_corr_metric_list   = cellfun(Original_Corr_Metric_fun, Tree_Metrics_cell(1, :));   

    % Correlation metric label
    if strcmp(correlation_metric, 'DBHqsm') == true
        correlation_metric_label = 'DBH';
    elseif strcmp(correlation_metric, 'UnmodDBHqsm') == true
        correlation_metric_label = 'unmod. DBH';
    else
        correlation_metric_label = regexprep(correlation_metric, '([A-Z])', ' ${lower($1)}');           % Inserts a space before each capital letter and makes them lower case

        if strcmp(correlation_metric_label(1), ' ')                                                     % This may result in a space being placed at the start, which is removed
            correlation_metric_label(1) = [];
        end

        correlation_metric_label = strrep(correlation_metric_label, 'trunk', 'stem');                   % Trunk is replaced by stem
    end



    % Correlation to the tree metrics' relative errors
    tree_metrics_cell   = Tree_Metrics.Values.Properties.VariableNames;
    number_tree_metrics = length(tree_metrics_cell);

    TreeQSM_corr_coeff_matrix   = zeros(number_beam_divergences, number_tree_metrics);
    TreeQSM_P_matrix            = zeros(number_beam_divergences, number_tree_metrics);
    FuzzyQSM_corr_coeff_matrix  = zeros(number_beam_divergences, number_tree_metrics);
    FuzzyQSM_P_matrix           = zeros(number_beam_divergences, number_tree_metrics);

    for m = 1 : number_tree_metrics
        % The analysed tree metric
        tree_metric = tree_metrics_cell{m};

        % Labels        
        if strcmp(tree_metric, 'DBHqsm')
            tree_metric_label = 'DBH';
        elseif strcmp(tree_metric, 'UnmodDBHqsm')
            tree_metric_label = 'unmod. DBH';
        else
            tree_metric_label = regexprep(tree_metric, '([A-Z])', ' ${lower($1)}');         % Inserts a space before each capital letter and makes them lower case
    
            if strcmp(tree_metric_label(1), ' ')                                            % This may result in a space being placed at the start, which is removed
                tree_metric_label(1) = [];
            end
    
            tree_metric_label = strrep(tree_metric_label, 'trunk', 'stem');                 % Trunk is replaced by stem
        end

        % The relative errors
        TreeQSM_Rel_Error_fun           = @(Tree_Metrics) abs(str2double(Tree_Metrics.Rel_Error('Discrete', tree_metric).Variables));
        TreeQSM_abs_rel_error_matrix    = cellfun(TreeQSM_Rel_Error_fun, Tree_Metrics_cell);

        FuzzyQSM_Rel_Error_fun          = @(Tree_Metrics) abs(str2double(Tree_Metrics.Rel_Error('Fuzzy', tree_metric).Variables));
        FuzzyQSM_abs_rel_error_matrix   = cellfun(FuzzyQSM_Rel_Error_fun, Tree_Metrics_cell);

        % Correlation coefficients
        for b = 1 : number_beam_divergences
            TreeQSM_abs_rel_error_list          = TreeQSM_abs_rel_error_matrix(b, :);
            [TreeQSM_corr_coeff, TreeQSM_P] = corrcoef(original_corr_metric_list, TreeQSM_abs_rel_error_list);
            TreeQSM_corr_coeff_matrix(b, m) = TreeQSM_corr_coeff(2);
            TreeQSM_P_matrix(b, m)          = TreeQSM_P(2);
    
            FuzzyQSM_abs_rel_error_list             = FuzzyQSM_abs_rel_error_matrix(b, :);
            [FuzzyQSM_corr_coeff, FuzzyQSM_P]   = corrcoef(original_corr_metric_list, FuzzyQSM_abs_rel_error_list);
            FuzzyQSM_corr_coeff_matrix(b, m)    = FuzzyQSM_corr_coeff(2);
            FuzzyQSM_P_matrix(b, m)             = FuzzyQSM_P(2);
        end

        % Plot
        figure_number = m;
        Fig = figure(figure_number);
        
        figure_name = sprintf('Correlation_Abs_Rel_Error_%s_%s', correlation_metric, tree_metric);
    
        set(Fig, 'name', figure_name, 'NumberTitle', 'off');
    
        % Set the size and white background color
        set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
        set(gcf, 'color', [1, 1, 1])    
        
        hold on
        grid on

        for b = 1 : number_beam_divergences
            % This beam divergence's data
            beam_divergence         = beam_divergence_list(b);
            beam_divergence_colour  = beam_divergence_colours{b};

            TreeQSM_abs_rel_error_list  = TreeQSM_abs_rel_error_matrix(b, :);
            FuzzyQSM_abs_rel_error_list = FuzzyQSM_abs_rel_error_matrix(b, :);

            scatter(original_corr_metric_list, TreeQSM_abs_rel_error_list, 'LineWidth', 2, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', beam_divergence_colour, 'DisplayName', sprintf('TreeQSM, %s = %.3g mrad', '\lambda', beam_divergence));
            scatter(original_corr_metric_list, FuzzyQSM_abs_rel_error_list, 'MarkerFaceColor', beam_divergence_colour, 'MarkerEdgeColor', beam_divergence_colour, 'DisplayName', sprintf('FuzzyQSM, %s = %.3g mrad', '\lambda', beam_divergence));
        end

        % Axes
        xlabel(sprintf('original %s [m]', correlation_metric_label));
        ylabel(sprintf('%s abs. rel. error [%%]', tree_metric_label));

        % Formatting
        legend('show', 'location', 'eastoutside');
        
        set(gca, 'FontSize', 20);
        set(gca, 'LineWidth', 2);
        
        hold off
        
        % Saving the figure    
        export_fig(figure_number, [figure_name, '.png']);
        movefile(sprintf('%s.png', figure_name), output_folder);
    end

%% Table %%
    % Method indices for consistency
    TreeQSM_ind     = 1;
    FuzzyQSM_ind    = 2;
    number_methods  = FuzzyQSM_ind;

    % Correlation coefficients and P-values for both methods and beam divergences
    corr_coeff_table_matrix = zeros(number_methods * number_beam_divergences, number_tree_metrics);  
    P_value_table_matrix    = zeros(number_methods * number_beam_divergences, number_tree_metrics);  
    row_labels              = cell(number_methods * number_beam_divergences, 1);

    for m = 1 : number_methods
        for b = 1 : number_beam_divergences
            row_ind         = sub2ind([number_methods, number_beam_divergences], m, b);

            if m == TreeQSM_ind
                corr_coeff_list = TreeQSM_corr_coeff_matrix(b, :);
                P_value_list    = TreeQSM_P_matrix(b, :);
                method_label    = 'TreeQSM';
            elseif m == FuzzyQSM_ind
                corr_coeff_list = FuzzyQSM_corr_coeff_matrix(b, :);
                P_value_list    = FuzzyQSM_P_matrix(b, :);
                method_label    = 'FuzzyQSM';
            end

            corr_coeff_table_matrix(row_ind, :) = corr_coeff_list;
            P_value_table_matrix(row_ind, :)    = P_value_list;

            beam_divergence     = beam_divergence_list(b);
            row_labels{row_ind} = sprintf('%s, lambda = %.3g mrad', method_label, beam_divergence);
        end
    end

    % Correlation coefficients table
    table_file_name = 'Corr_Coeff_Abs_Rel_Error.xls';
    Table_Formatter(corr_coeff_table_matrix, num_decimal_digits, output_format, row_labels, tree_metrics_cell, table_file_name, Print);
    movefile(table_file_name, output_folder);

    % P-values table
    table_file_name = 'P_Values_Abs_Rel_Error.xls';
    Table_Formatter(P_value_table_matrix, num_decimal_digits, output_format, row_labels, tree_metrics_cell, table_file_name, Print);
    movefile(table_file_name, output_folder);
