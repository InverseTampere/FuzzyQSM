% If the uncertainty of TreeQSM with and without resampling has been determined, this script subtracts them to estimate the TreeQSM stochasticity

clear variables
close all
clc

%% Inputs %%
    % Folder names
    resampling_folder       = 'Small_Tree_Results';
    no_resampling_folder    = 'Small_Tree_Results_No_Resampling';

    % Table settings
    num_decimal_digits          = 3;            % [-]
    output_format               = 'Float';      % [Float, Integer, Exponential]

%% Retrieve the data %%
    % Height data
    QSM_Height_Data_Res         = load(sprintf('%s/%s', resampling_folder, 'QSM_Height_Data.mat'));
    QSM_Height_Data_NonRes      = load(sprintf('%s/%s', no_resampling_folder, 'QSM_Height_Data.mat'));
    scanner_combinations_cell   = QSM_Height_Data_Res.scanner_combinations_cell;
    number_scanner_combinations = length(scanner_combinations_cell);

%% Save folder %%
    % A folder is created with the current time as its name to save the results in
    current_time        = datetime('now', 'format', 'yyyyMMdd_HHmmss');                         % Date followed by clock time       
    save_folder_name    = string(current_time);                                                 % Converted to string array that acts as the name of the folder

    mkdir(save_folder_name);

    % String labels for the scanner combinations and their colour map
    scanner_combination_labels  = cellfun(@num2str, scanner_combinations_cell, 'UniformOutput', false);
    scanner_comb_cmap           = cbrewer('qual', 'Set1', max(number_scanner_combinations, 3));

%% Tree parameters %%
    for c = 1 : number_scanner_combinations
        % Data 
        scanner_combination = scanner_combinations_cell{c};
        table_file_name     = sprintf('Tree_Parameter_%s_Results.xls', num2str(scanner_combination));
        
        Res_Table       = readtable(sprintf('%s/%s', resampling_folder, table_file_name), 'VariableNamingRule', 'preserve');
        Non_Res_Table   = readtable(sprintf('%s/%s', no_resampling_folder, table_file_name), 'VariableNamingRule', 'preserve');

        % Difference in mean values and standard deviations
        res_mean_cell       = Res_Table.('MC mean');
        res_mean_list       = cellfun(@str2double, res_mean_cell);
        non_res_mean_cell   = Non_Res_Table.('MC mean');
        non_res_mean_list   = cellfun(@str2double, non_res_mean_cell);

        mean_delta_list = res_mean_list - non_res_mean_list;
        
        res_std_cell        = Res_Table.('MC std');
        res_std_list        = cellfun(@str2double, res_std_cell);
        non_res_std_cell    = Non_Res_Table.('MC std');
        non_res_std_list    = cellfun(@str2double, non_res_std_cell);

        std_delta_list = res_std_list - non_res_std_list;

        % Table
        tree_parameter_fields   = Res_Table.Row;        
        column_names            = {'mean res', 'mean (res - non-res)', 'std res', 'std (res - non-res)'};
        tree_table_data         = [res_mean_list, mean_delta_list, res_std_list, std_delta_list];
        Print                   = true;
        Table_Formatter(tree_table_data, num_decimal_digits, output_format, tree_parameter_fields, column_names, table_file_name, Print);

        movefile(table_file_name, save_folder_name);
    end

%% Height-based results %%
    % Data
    bin_edge_heights    = QSM_Height_Data_Res.bin_edge_heights;
    Height_Data_Res     = QSM_Height_Data_Res.QSM_Height_Data;
    Height_Data_NonRes  = QSM_Height_Data_NonRes.QSM_Height_Data;

    % Average height for each bin
    bin_avg_heights = (bin_edge_heights(2 : end) + bin_edge_heights(1 : end - 1)) / 2;

    % Plotted data
    height_metrics_cell         = {'branch_volume', 'branch_number_cylinders', 'stem_volume', 'stem_number_cylinders'};
    height_metric_units_cell    = {'m^3', '-', 'm^3', '-'};
    number_metrics              = length(height_metrics_cell); 

    for m = 1 : number_metrics
        % The plotted height metric
        height_metric   = height_metrics_cell{m};
        metric_unit     = height_metric_units_cell{m};

        % Data types
        data_fields         = fieldnames(Height_Data_Res(1).(height_metric));
        number_data_fields  = length(data_fields);

        for d = 1 : number_data_fields
            % Plotting the data
            data_field = data_fields{d};

            figure_number = 2000 + 100*m + d*10;
            figure(figure_number)

            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.7 0.9])
            set(gcf, 'color', [1, 1, 1])   

            hold on
            grid on

            % Each scanner combination's data
            for c = 1 : number_scanner_combinations
                % The data
                scanner_comb_colour = scanner_comb_cmap(c, :);
                scanner_combination = scanner_combinations_cell{c};
                resampled_data      = Height_Data_Res(c).(height_metric).(data_field);
                non_resampled_data  = Height_Data_NonRes(c).(height_metric).(data_field);

                % Plot
                plot(resampled_data - non_resampled_data, bin_avg_heights, 'LineWidth', 2, 'Color', scanner_comb_colour, 'DisplayName', num2str(scanner_combination));
            end

            % Axes
            metric_label = strrep(height_metric, '_', ' ');
            
            if strcmp(data_field, 'mean')
                field_label = '\mu';
                data_unit   = metric_unit;
            elseif strcmp(data_field, 'std')
                field_label = '\sigma';
                data_unit   = metric_unit;
            elseif strcmp(data_field, 'relative_uncertainty')
                field_label = 'rel. \sigma';
                data_unit   = '-';
            end

            x_axis_label = sprintf('%s R - NR %s [%s]', metric_label, field_label, data_unit);
            xlabel(x_axis_label);

            ylabel('H [m]');

            % Legend and text
            legend('show', 'location', 'eastoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    

            % Saving the figure
            figure_name = sprintf('TreeQSM_Height_Sampling_%s_%s', height_metric, data_field);
            export_fig(figure_number, [figure_name, '.fig']);
            export_fig(figure_number, [figure_name, '.png']);

            movefile([figure_name, '*'], save_folder_name); 
        end
    end