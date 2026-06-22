% This script plots the fitted versus true tree metric values for all trees and both beam divergences
% Note that the results need to already have been generated

clear variables
close all
beep off
clc

%% Inputs %%
    % Beam divergence folders
    beam_divergence_1_folder    = 'Data/Results/030_beamdivergence';
    beam_divergence_2_folder    = 'Data/Results/0175_beamdivergence';

    % Output folder
    output_folder               = 'Data/Results/Tree_Metrics';

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

    % Tree metrics tables
    Beam_Divergence_Tree_Metrics_cell = cell(1, number_beam_divergences);

    for b = 1 : number_beam_divergences
        % This beam divergence's folders
        beam_divergence_folder = beam_divergence_folders{b};

        tree_folders    = dir(beam_divergence_folder);
        tree_folders    = {tree_folders(:).name};
    
        dot_folders     = contains(tree_folders, '.');                        % Remove the . and .. folders
        tree_folders    = tree_folders(~dot_folders);
    
        rel_uncert_bool = strcmp(tree_folders, 'Relative_Uncertainty');       % Remove the relative uncertainty results if they exist
        tree_folders    = tree_folders(~rel_uncert_bool);
    
        stem_fit_bool   = strcmp(tree_folders, 'Stem_Fitting');               % Remove the stem fitting results if they exist
        tree_folders    = tree_folders(~stem_fit_bool);        
        number_trees    = length(tree_folders);

        % The tree metrics
        Tree_Metrics_cell = cell(1, number_trees);

        for t = 1 : number_trees
            % This tree's data
            tree_folder = tree_folders{t};

            results_file_name       = sprintf('%s/%s/FuzzyQSM_Outputs.mat', beam_divergence_folder, tree_folder);
            Results_File            = load(results_file_name);
            Tree_Metrics            = Results_File.Results_Tables.Tree_Metrics;
            Tree_Metrics_cell{t}    = Tree_Metrics;
        end

        Beam_Divergence_Tree_Metrics_cell{b} = Tree_Metrics_cell;
    end

    % Determine the species of each of the trees
    Tree_Folder_Split_fun   = @(tree_folder) strsplit(tree_folder, '_');
    tree_folder_parts_cell  = cellfun(Tree_Folder_Split_fun, tree_folders, 'UniformOutput', false);
    Tree_Species_fun        = @(tree_folder_parts) tree_folder_parts{1};
    tree_species_cell       = cellfun(Tree_Species_fun, tree_folder_parts_cell, 'UniformOutput', false);

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

%% Plots %%
    % Colour map for each species
    species_cmap    = cbrewer('qual', 'Set1', max(3, number_species));

    % Marker for each beam divergence
    marker_cell     = {'o', '^', '+', 'pentagram'};

    % Tree metrics data
    QSM_labels                              = {'Discrete', 'Fuzzy', 'Original'};
    [Discr_ind, Fuzzy_ind, Original_ind]    = Column_Deal(1: 3);
    number_QSMs                             = length(QSM_labels);
    tree_metrics_cell                       = cell(number_beam_divergences, number_QSMs);

    for b = 1 : number_beam_divergences
        % This beam divergence's data
        Tree_Metrics_cell = Beam_Divergence_Tree_Metrics_cell{b};

        for q = 1 : number_QSMs
            % This QSM's tree metrics data
            QSM_label               = QSM_labels{q};
            Tree_Metrics_fun        = @(Tree_Metrics) cellfun(@str2double, Tree_Metrics.Values{QSM_label, :});
            tree_metrics_cell_q     = cellfun(Tree_Metrics_fun, Tree_Metrics_cell, 'UniformOutput', false);
            tree_metrics_matrix_q   = vertcat(tree_metrics_cell_q{:});
            tree_metrics_cell{b, q} = tree_metrics_matrix_q;
        end
    end

    total_metric_matrix = vertcat(tree_metrics_cell{:});

    %--% Plots for each tree metric %--%
    tree_metrics        = Tree_Metrics.Rel_Error.Properties.VariableNames;
    number_tree_metrics = length(tree_metrics);

    for t = 1 : number_tree_metrics
        % This metric's data
        tree_metric         = tree_metrics{t};
        total_metric_list   = total_metric_matrix(:, t);
        
        % Data bounds with margin
        data_LB     = min(total_metric_list);
        max_value   = max(total_metric_list);
        data_ampl   = max_value - data_LB;
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
            % DBH or height
            tree_metric_unit = 'm';

            if strcmp(tree_metric, 'DBHqsm')
                tree_metric_label = 'DBH';
            elseif strcmp(tree_metric, 'UnmodDBHqsm')
                tree_metric_label = 'unmod. DBH';
            elseif strcmp(tree_metric, 'TreeHeight')
                tree_metric_label = 'tree height';
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

        % Fuzzy and discrete legend entries for the two beam divergences
        for b = 1 : number_beam_divergences
            beam_divergence = beam_divergence_list(b);
            marker          = marker_cell{b};

            scatter(NaN, NaN, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'marker', marker, 'DisplayName', sprintf('TreeQSM, %s = %.3f mrad', '\lambda', beam_divergence));
            scatter(NaN, NaN, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'marker', marker, 'DisplayName', sprintf('FuzzyQSM, %s = %.3f mrad', '\lambda', beam_divergence));
        end

        % 0% error level
        plot([data_LB, data_UB], [data_LB, data_UB], 'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'Off');     

        if contains(tree_metric, 'Volume')
            % 20% error level
            plot([data_LB, data_UB], 0.80*[data_LB, data_UB], 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':', 'HandleVisibility', 'Off');        
            plot([data_LB, data_UB], 1.20*[data_LB, data_UB], 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':', 'HandleVisibility', 'Off');       
        else
            % 10% error level
            plot([data_LB, data_UB], 0.90*[data_LB, data_UB], 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':', 'HandleVisibility', 'Off');        
            plot([data_LB, data_UB], 1.10*[data_LB, data_UB], 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':', 'HandleVisibility', 'Off');       
        end

        % Fuzzy and discrete values vs. the original values
        for s = 1 : number_species
            % This species' data
            species_colour  = species_cmap(s, :);
            species_bool    = species_bool_cell{s};
            tree_species    = unique_tree_species{s};

            for b = 1 : number_beam_divergences
                % This beam divergence's data
                marker                      = marker_cell{b};
                original_tree_metrics_list  = tree_metrics_cell{b, Original_ind}(species_bool, t);
                discrete_tree_metrics_list  = tree_metrics_cell{b, Discr_ind}(species_bool, t);
                fuzzy_tree_metrics_list     = tree_metrics_cell{b, Fuzzy_ind}(species_bool, t);

                % Scatter plots
                scatter(original_tree_metrics_list, discrete_tree_metrics_list, 'marker', marker, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', species_colour, 'LineWidth', 2, 'HandleVisibility', 'Off');
                scatter(original_tree_metrics_list, fuzzy_tree_metrics_list, 'marker', marker, 'MarkerFaceColor', species_colour, 'MarkerEdgeColor', species_colour, 'LineWidth', 2, 'HandleVisibility', 'Off');
            end
        end

        % Axes
        if contains(tree_metric, 'Volume')
            % Switch to log-log plot
            set(gca, 'XScale', 'log', 'YScale', 'log')      
        end

        xlim([data_LB, data_UB]);
        ylim([data_LB, data_UB]);

        xlabel(sprintf('true %s [%s]', tree_metric_label, tree_metric_unit));
        ylabel(sprintf('fitted %s [%s]', tree_metric_label, tree_metric_unit));

        % Legend
        legend('show', 'location', 'eastoutside');

        set(gca, 'FontSize', 20);
        set(gca, 'LineWidth', 2);

        hold off

        % Saving the figure        
        export_fig(figure_number, sprintf('%s.fig', figure_name));
        export_fig(figure_number, sprintf('%s.png', figure_name), '-transparent');

        movefile([figure_name, '*'], output_folder); 
    end

    %--% Separate species legend %--%
    figure_number   = figure_number + 1;
    Fig             = figure(figure_number);
    figure_name     = 'Species_Legend';

    set(Fig, 'name', figure_name, 'NumberTitle', 'off');

    % Set the size and white background color
    set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
    set(gcf, 'color', [1, 1, 1])    

    hold on

    % Species entries
    for s = 1 : number_species
        species_colour  = species_cmap(s, :);
        tree_species    = unique_tree_species{s};

        scatter(NaN, NaN, 'MarkerEdgeColor', species_colour, 'MarkerFaceColor', species_colour, 'DisplayName', tree_species);
    end

    % Legend
    legend('show', 'location', 'north', 'orientation', 'horizontal');
    set(gca, 'FontSize', 15);
    set(gca, 'LineWidth', 2);

    axis off

    export_fig(figure_number, sprintf('%s.fig', figure_name));
    export_fig(figure_number, sprintf('%s.png', figure_name), '-transparent');

    movefile([figure_name, '*'], output_folder); 

