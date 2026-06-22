% This script creates plots of the relative uncertainty results for the two beam divergences
% The results themselves should have been generated already

clear variables
close all
clc

%% Inputs %%
    % Beam divergence folders
    beam_divergence_1_folder    = 'Data/Results/030_beamdivergence/Relative_Uncertainty';
    beam_divergence_2_folder    = 'Data/Results/0175_beamdivergence/Relative_Uncertainty';

    % Output folder
    output_folder               = 'Data/Results/Relative_Uncertainty';

    % Plot settings
    axes_width                  = 0.48;      % Width to ensure consistent sizing regardless of legend size

%% Retrieve the data %%
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

    % Data
    [Tree_Metrics_cell, fuzzy_cyl_percentage_cell, fuzzy_points_percentage_cell] = deal(cell(1, number_beam_divergences));

    for b = 1 : number_beam_divergences
        % This beam divergence's data
        beam_divergence_folder  = beam_divergence_folders{b};
        results_file_path       = sprintf('%s/Rel_Uncertainty_Results.mat', beam_divergence_folder);
        Results_File            = load(results_file_path);
        
        Results_Tables_cell             = Results_File.Results_Tables_cell;
        fuzzy_cyl_percentage_cell{b}    = Results_File.fuzzy_cyl_percentage_matrix;
        fuzzy_points_percentage_cell{b} = Results_File.fuzzy_points_percentage_matrix;
        rel_uncertainty_threshold_steps = Results_File.rel_uncertainty_threshold_steps;

        max_threshold   = max(rel_uncertainty_threshold_steps);
        number_steps    = length(rel_uncertainty_threshold_steps);

        % The tree metrics
        Tree_Metrics_fun        = @(Results_Tables) Results_Tables.Tree_Metrics;
        Tree_Metrics_cell{b}    = cellfun(Tree_Metrics_fun, Results_Tables_cell, 'UniformOutput', false);
    end

%% Diagnostic metrics %%
    % Method colours
    set_2_cmap      = cbrewer('qual', 'Set2', 3);
    TreeQSM_colour  = set_2_cmap(1, :);
    set_1_cmap      = cbrewer('qual', 'Set1', 3);
    FuzzyQSM_colour = set_1_cmap(1, :);

    beam_divergence_linestyles = {'-', '--'};

    %--% Percentage of fuzzy cylinders %--%
    figure_number   = 1;
    Fig             = figure(figure_number);
    figure_name     = 'Fuzzy_Cylinders_Percentage';

    set(Fig, 'name', figure_name, 'NumberTitle', 'off');

    % Set the size and white background color
    set(gcf, 'Units', 'Normalized', 'Position', [0 0 1 1]);
    set(gcf, 'color', [1, 1, 1]);
    
    hold on
    grid on

    for b = 1 : number_beam_divergences
        beam_divergence_linestyle       = beam_divergence_linestyles{b};
        beam_divergence                 = beam_divergence_list(b);
        fuzzy_cyl_percentage_matrix     = fuzzy_cyl_percentage_cell{b};
        mean_fuzzy_cyl_percentage_list  = mean(fuzzy_cyl_percentage_matrix, 1);
        
        plot(rel_uncertainty_threshold_steps, mean_fuzzy_cyl_percentage_list, 'LineWidth', 3, 'color', 'b', 'LineStyle', beam_divergence_linestyle, 'DisplayName', sprintf('mean, %s = %.3f mrad', '\lambda', beam_divergence));
    end

    % Axes
    xlim([0, max_threshold]);
    xlabel(sprintf('%s [-]', '\nu'));
    ylabel('fuzzy fitted cylinders [%]');

    % Formatting
    legend('show', 'location', 'eastoutside');
    set(gca, 'FontSize', 25);
    set(gca, 'LineWidth', 3);

    axes_position       = get(gca, 'Position');
    axes_position(3)    = axes_width;
    set(gca, 'Position', axes_position)

    hold off

    % Saving the figure     
    export_fig(figure_number, sprintf('%s.fig', figure_name));
    export_fig(figure_number, sprintf('%s.png', figure_name));
    movefile([figure_name, '*'], output_folder); 

    %--% Percentage of fuzzy points %--%
    figure_number   = figure_number + 1;
    Fig             = figure(figure_number);
    figure_name     = 'Fuzzy_Points_Percentage';

    set(Fig, 'name', figure_name, 'NumberTitle', 'off');

    % Set the size and white background color
    set(gcf, 'Units', 'Normalized', 'Position', [0 0 1 1]);
    set(gcf, 'color', [1, 1, 1]);
    
    hold on
    grid on

    for b = 1 : number_beam_divergences
        beam_divergence_linestyle           = beam_divergence_linestyles{b};
        beam_divergence                     = beam_divergence_list(b);
        fuzzy_points_percentage_matrix      = fuzzy_points_percentage_cell{b};
        mean_fuzzy_points_percentage_list   = mean(fuzzy_points_percentage_matrix, 1);

        plot(rel_uncertainty_threshold_steps, mean_fuzzy_points_percentage_list, 'LineWidth', 3, 'color', 'b', 'LineStyle', beam_divergence_linestyle, 'DisplayName', sprintf('mean, %s = %.3f mrad', '\lambda', beam_divergence));
    end

    % Axes
    xlim([0, max_threshold]);
    xlabel(sprintf('%s [-]', '\nu'));
    ylabel('points for fuzzy fitting [%]');

    % Formatting
    legend('show', 'location', 'eastoutside');

    set(gca, 'FontSize', 25);
    set(gca, 'LineWidth', 3);

    axes_position       = get(gca, 'Position');
    axes_position(3)    = axes_width;
    set(gca, 'Position', axes_position)

    hold off

    % Saving the figure     
    export_fig(figure_number, sprintf('%s.fig', figure_name));
    export_fig(figure_number, sprintf('%s.png', figure_name));
    movefile([figure_name, '*'], output_folder); 

%% Tree metric errors %%
    % Analysed tree metrics
    tree_metrics        = Tree_Metrics_cell{1}{1}.Rel_Error.Properties.VariableNames;
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

        %--% Figure %--%
        figure_number   = figure_number + 1;
        Fig             = figure(figure_number);
        figure_name     = sprintf('Relative_Error_%s', tree_metric);
  
        set(Fig, 'name', figure_name, 'NumberTitle', 'off');

        % Set the size and white background color
        set(gcf, 'Units', 'Normalized', 'Position', [0 0 1 1]);
        set(gcf, 'color', [1, 1, 1]);

        hold on
        grid on
        
        for b = 1 : number_beam_divergences
            % This beam divergence's data
            beam_divergence_linestyle           = beam_divergence_linestyles{b};
            beam_divergence                     = beam_divergence_list(b);

            Discr_Tree_Metrics_fun              = @(Tree_Metrics) str2double(Tree_Metrics.Rel_Error('Discrete', tree_metric).Variables);
            discr_tree_metrics_rel_error_list   = cellfun(Discr_Tree_Metrics_fun, Tree_Metrics_cell{b}(:, 1), 'UniformOutput', true);
            discrete_mean_rel_error             = mean(discr_tree_metrics_rel_error_list);
    
            Fuzzy_Tree_Metrics_fun              = @(Tree_Metrics) str2double(Tree_Metrics.Rel_Error('Fuzzy', tree_metric).Variables);
            fuzzy_tree_metrics_rel_error_matrix = cellfun(Fuzzy_Tree_Metrics_fun, Tree_Metrics_cell{b}, 'UniformOutput', true);
            fuzzy_mean_rel_error_list           = mean(fuzzy_tree_metrics_rel_error_matrix, 1);

            % Discrete error
            plot(rel_uncertainty_threshold_steps, discrete_mean_rel_error * ones(1, number_steps), 'LineWidth', 3, 'LineStyle', beam_divergence_linestyle, 'color', TreeQSM_colour, 'DisplayName', sprintf('TreeQSM mean, %s = %.3f mrad', '\lambda', beam_divergence));
    
            % Fuzzy error
            plot(rel_uncertainty_threshold_steps, fuzzy_mean_rel_error_list, 'LineWidth', 3, 'LineStyle', beam_divergence_linestyle, 'color', FuzzyQSM_colour, 'DisplayName', sprintf('MixedQSM mean, %s = %.3f mrad', '\lambda', beam_divergence));
        end

        % Axes
        xlim([0, max_threshold]);
	    xlabel(sprintf('%s [-]', '\nu'));
        ylabel(sprintf('%s rel. error [%%]', tree_metric_label));

        % Formatting
        legend('show', 'location', 'eastoutside');
        set(gca, 'FontSize', 25);
        set(gca, 'LineWidth', 3);

        axes_position       = get(gca, 'Position');
        axes_position(3)    = axes_width;
        set(gca, 'Position', axes_position)

        hold off

        % Saving the figure     
        export_fig(figure_number, sprintf('%s.fig', figure_name));
        export_fig(figure_number, sprintf('%s.png', figure_name));
        movefile([figure_name, '*'], output_folder); 
    end
