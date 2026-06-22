% This script looks at the stem fitting results for two beam divergence angles and plots them together

clear variables
close all
clc

%% Inputs %%
    % Beam divergence folders
    beam_divergence_1_folder    = 'Data/Results/030_beamdivergence/Stem_Fitting';
    beam_divergence_2_folder    = 'Data/Results/0175_beamdivergence/Stem_Fitting';

    beam_divergence_linestyles  = {'-', '--'};
    
    % Output folder
    output_folder               = 'Data/Results/Stem_Fitting';

    % Plot settings
    axes_width                  = 0.48;

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

    % Results
    [Tree_Metrics_Rel_Errors, comp_time_cell, percent_distr_cell] = deal(cell(1, number_beam_divergences));

    for b = 1 : number_beam_divergences
        % This beam divergence's data
        beam_divergence_folder  = beam_divergence_folders{b};
        results_file_path       = sprintf('%s/Tree_Metrics_Rel_Errors.mat', beam_divergence_folder);

        Results_File                    = load(results_file_path);
        Tree_Metrics_Rel_Errors{b}      = Results_File.Tree_Metrics_Rel_Errors_cell;
        comp_time_cell{b}               = Results_File.comp_time_matrix;
        percent_distr_cell{b}           = Results_File.percent_distr_matrix;
        max_number_distributions_list   = Results_File.max_number_distributions_list;    
        number_steps                    = length(max_number_distributions_list);
    end

%% Diagnostic metrics %%
    %--% Computational time %--%
    figure_number = 1;
    Fig = figure(figure_number);
    
    figure_name = 'Num_Distr_Comp_Time';

    set(Fig, 'name', figure_name, 'NumberTitle', 'off');

    % Set the size and white background color
    set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 1, 0.9])
    set(gcf, 'color', [1, 1, 1])    
    
    hold on
    grid on

    for b = 1 : number_beam_divergences
        % This beam divergence's relative computational time
        beam_divergence             = beam_divergence_list(b);
        beam_divergence_linestyle   = beam_divergence_linestyles{b};
        comp_time_matrix            = comp_time_cell{b};

        rel_comp_time_matrix    = comp_time_matrix ./ comp_time_matrix(:, 1) * 100;
        mean_rel_comp_time_list = mean(rel_comp_time_matrix, 1);

        plot(1:2, mean_rel_comp_time_list(1:2), 'LineWidth', 1, 'LineStyle', ':', 'Color', 'b', 'HandleVisibility', 'Off');
        plot(2:number_steps, mean_rel_comp_time_list(2:number_steps), 'LineWidth', 3, 'LineStyle', beam_divergence_linestyle, 'Color', 'b', 'DisplayName', sprintf('mean, %s = %.3f mrad', '\lambda', beam_divergence));
    end

    % Axes
    if number_steps ~= 1
        xlim([1, number_steps]);
        xticks(1 : number_steps);
        xlabel('max. nr. distributions [-]');
        xticklabels(['All', string(max_number_distributions_list(2 : number_steps))]);
    else
        xlabel('max nr. distributions');
    end
    
    ylabel('rel. computational time [%]');

    % Formatting
    legend('show', 'location', 'eastoutside');

    set(gca, 'FontSize', 25);
    set(gca, 'LineWidth', 3);

    axes_position       = get(gca, 'Position');
    axes_position(3)    = axes_width;
    set(gca, 'Position', axes_position)
    
    hold off
    
    % Saving the figure    
    export_fig(figure_number, [figure_name, '.png']);
    movefile(sprintf('%s.png', figure_name), output_folder);

    %--% Percentage of remaining points/distributions %--%
    figure_number = 2;
    Fig = figure(figure_number);
    
    figure_name = 'Num_Distr_Percent_Distr';

    set(Fig, 'name', figure_name, 'NumberTitle', 'off');

    % Set the size and white background color
    set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 1, 0.9])
    set(gcf, 'color', [1, 1, 1])    
    
    hold on
    grid on

    for b = 1 : number_beam_divergences
        % This beam divergence's number of fuzzy points
        beam_divergence             = beam_divergence_list(b);
        beam_divergence_linestyle   = beam_divergence_linestyles{b};
        percent_distr_matrix        = percent_distr_cell{b};
        mean_percent_distr_list     = mean(percent_distr_matrix, 1);

        plot(1:2, mean_percent_distr_list(1:2), 'LineWidth', 1, 'Color', 'b', 'LineStyle', ':', 'HandleVisibility', 'Off');
        plot(2:number_steps, mean_percent_distr_list(2:number_steps), 'LineWidth', 3, 'Color', 'b', 'LineStyle', beam_divergence_linestyle, 'DisplayName', sprintf('mean, %s = %.3f mrad', '\lambda', beam_divergence));
    end

    % Axes
    if number_steps ~= 1
        xlim([1, number_steps]);
        xticks(1 : number_steps);
        xlabel('max. nr. distributions [-]');
        xticklabels(['All', string(max_number_distributions_list(2 : number_steps))]);
    else
        xlabel('max nr. distributions');
    end
    
    ylabel('remaining points [%]');

    % Formatting
    legend('show', 'location', 'eastoutside');

    set(gca, 'FontSize', 25);
    set(gca, 'LineWidth', 3);

    axes_position       = get(gca, 'Position');
    axes_position(3)    = axes_width;
    set(gca, 'Position', axes_position)

    
    hold off
    
    % Saving the figure    
    export_fig(figure_number, [figure_name, '.png']);
    movefile(sprintf('%s.png', figure_name), output_folder);

%% Tree metric errors %%
    % Metrics
    tree_metrics    = {'TrunkVolume', 'DBHqsm'};           % The others don't make sense to analyse as no branches were fitted
    metric_labels   = {'stem volume', 'DBH'};
    number_metrics  = length(tree_metrics);

    % Method colours
    set_2_cmap      = cbrewer('qual', 'Set2', 3);
    TreeQSM_colour  = set_2_cmap(1, :);
    set_1_cmap      = cbrewer('qual', 'Set1', 3);
    FuzzyQSM_colour = set_1_cmap(1, :);

    for m = 1 : number_metrics
        % This metric's data
        tree_metric     = tree_metrics{m};
        metric_label    = metric_labels{m};
        
        % Plot
        figure_number = figure_number + 1;
        Fig = figure(figure_number);
        
        figure_name = sprintf('%s_Num_Distr_Rel_Errors', tree_metric);

        set(Fig, 'name', figure_name, 'NumberTitle', 'off');
    
        % Set the size and white background color
        set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 1, 0.9])
        set(gcf, 'color', [1, 1, 1])    
        
        hold on
        grid on
        
        % Zero line for clarity
        plot([1, number_steps], [0, 0], 'LineWidth', 3, 'color', 'k', 'DisplayName', 'zero-line', 'HandleVisibility', 'Off');

        for b = 1 : number_beam_divergences
            % This beam divergence's data
            beam_divergence                 = beam_divergence_list(b);
            beam_divergence_linestyle       = beam_divergence_linestyles{b};
            Tree_Metrics_Rel_Errors_cell    = Tree_Metrics_Rel_Errors{b};

            % Discrete relative errors for each tree
            Discr_Rel_Error_fun     = @(Tree_Metrics_Rel_Error) str2double(Tree_Metrics_Rel_Error('Discrete', tree_metric).Variables);
            discr_rel_error_list    = cellfun(Discr_Rel_Error_fun, Tree_Metrics_Rel_Errors_cell(:, 1));                     
            discr_mean              = mean(discr_rel_error_list);
    
            % Fuzzy relative errors for each tree and max. number of distributions
            Fuzzy_Rel_Error_fun     = @(Tree_Metrics_Rel_Error) str2double(Tree_Metrics_Rel_Error('Fuzzy', tree_metric).Variables);
            fuzzy_rel_error_matrix  = cellfun(Fuzzy_Rel_Error_fun, Tree_Metrics_Rel_Errors_cell);
            fuzzy_mean_list         = mean(fuzzy_rel_error_matrix, 1);

            % Discrete mean rel. error 
            plot([1, number_steps], discr_mean * [1, 1], 'LineWidth', 3, 'LineStyle', beam_divergence_linestyle, 'Color', TreeQSM_colour, 'DisplayName', sprintf('TreeQSM mean, %s = %.3f mrad', '\lambda', beam_divergence));
    
            % Fuzzy mean rel. error
            plot(1 : 2, fuzzy_mean_list(1:2), 'LineWidth', 1, 'Color', FuzzyQSM_colour, 'LineStyle', ':', 'HandleVisibility', 'Off');
            plot(2 : number_steps, fuzzy_mean_list(2:end), 'LineWidth', 3, 'Color', FuzzyQSM_colour, 'LineStyle', beam_divergence_linestyle, 'DisplayName', sprintf('FuzzyQSM mean, %s = %.3f mrad', '\lambda', beam_divergence));
        end

        % Axes
        if number_steps ~= 1
            xlim([1, number_steps]);
            xticks(1 : number_steps);
            xlabel('max. nr. distributions [-]');
            xticklabels(['All', string(max_number_distributions_list(2 : number_steps))]);
        else
            xlabel('max nr. distributions');
        end
    
        ylabel(sprintf('%s rel. error [%%]', metric_label));
    
        % Legend
        legend('show', 'location', 'eastoutside');

        % Formatting
        set(gca, 'FontSize', 25);
        set(gca, 'LineWidth', 3);

        axes_position       = get(gca, 'Position');
        axes_position(3)    = axes_width;
        set(gca, 'Position', axes_position)

        hold off
        
        % Saving the figure    
        export_fig(figure_number, [figure_name, '.png']);
        movefile(sprintf('%s.png', figure_name), output_folder);
    end