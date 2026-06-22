% This script compares the results fo two sets of QSMs, i.e. the mean and relative uncertainty of various metrics
% Note that the relative difference is 2*(val_a - val_b) / (val_a + val_b)

function TreeUQSM_Comparison(MC_QSM_cell_a, MC_QSM_cell_b, QSM_labels_cell, scanner_combinations_cell)

    %% Manual i nputs %%
        % Height binning
        height_step                 = 0.2;                                      % [m] For the height-based results
    
        % Plot settings
        Plot                        = false;                                     % [true, false] Whether or not figures are plotted
        Separate_Stem               = true;                                     % [true, false] Separate the stem from the branches for the plots
    
        % Table settings
        Print                       = true;                                     % [true, false] Whether or not tables are displayed
        num_decimal_digits          = 5;                                        % [-] Number of decimal digits of the table data
        output_format               = 'Float';                                  % [Float, Integer, Exponential] Format of the table data

    %% Preliminary %%
        % A folder is created with the current time as its name to save the results in
        current_time        = datetime('now', 'format', 'yyyyMMdd_HHmmss');                         % Date followed by clock time       
        save_folder_name    = string(current_time);                                                 % Converted to string array that acts as the name of the folder
    
        mkdir(save_folder_name);

        % Scanner combination colour map
        scanner_combination_labels  = cellfun(@num2str, scanner_combinations_cell, 'UniformOutput', false);
        number_scanner_combinations = length(scanner_combinations_cell);
        scanner_combination_cmap    = cbrewer('qual', 'Set1', max(number_scanner_combinations, 3));

    %% Remove non-existant branches %%
        % Function handle to trim the QSMs
        Trimming_fun = @(QSM) TreeQSM_Branch_Trimmer(QSM);

        % Monte Carlo QSMs
        for c = 1 : number_scanner_combinations
            % First set of QSMs
            QSM_cell_a = MC_QSM_cell_a{c};
            QSM_cell_a = cellfun(Trimming_fun, QSM_cell_a, 'UniformOutput', false);

            MC_QSM_cell_a{c} = QSM_cell_a;

            % Second set of QSMs
            QSM_cell_b = MC_QSM_cell_b{c};
            QSM_cell_b = cellfun(Trimming_fun, QSM_cell_b, 'UniformOutput', false);

            MC_QSM_cell_b{c} = QSM_cell_b;
        end

    %% Relative volume %%
        % Function handle for the relative volume of each branch
        Relative_Volume_fun = @(QSM) TreeQSM_Relative_Volume(QSM);

        % Monte Carlo QSMs
        for c = 1 : number_scanner_combinations
            % First set of QSMs
            QSM_cell_a = MC_QSM_cell_a{c};
            QSM_cell_a = cellfun(Relative_Volume_fun, QSM_cell_a, 'UniformOutput', false);

            MC_QSM_cell_a{c} = QSM_cell_a;

            % Second set of QSMs
            QSM_cell_b = MC_QSM_cell_b{c};
            QSM_cell_b = cellfun(Relative_Volume_fun, QSM_cell_b, 'UniformOutput', false);

            MC_QSM_cell_b{c} = QSM_cell_b;
        end

    %% Branch parameters %%
        % Analysed parameters
        branch_parameter_fields     = {'diameter', 'length', 'volume', 'RelativeVolume'};
        branch_parameter_units      = {'m', 'm', 'm^3', '%'};
        number_branch_parameters    = length(branch_parameter_fields);
    
        % Their data
        [Branch_Order_Parameters_a, Number_Branches_a, max_branch_order_a]  = TreeUQSM_Branch_Order_Parameters([], MC_QSM_cell_a, branch_parameter_fields);
        [Branch_Order_Parameters_b, Number_Branches_b, max_branch_order_b]  = TreeUQSM_Branch_Order_Parameters([], MC_QSM_cell_b, branch_parameter_fields);
        max_branch_order                                                    = min(max_branch_order_a, max_branch_order_b);
    
        if Separate_Stem == true
            branch_order_cell = {0, 1:max_branch_order};
        else
            branch_order_cell = {0:max_branch_order};
        end
    
        number_branch_sets = length(branch_order_cell);
    
        table_row_labels = string(0 : max_branch_order);
    
        % Results for each parameter
        for p = 1 : number_branch_parameters
            % This parameter's data
            branch_parameter        = branch_parameter_fields{p};
            branch_parameter_unit   = branch_parameter_units{p};
            Branch_Parameter_Data_a = Branch_Order_Parameters_a.(branch_parameter);
            Branch_Parameter_Data_b = Branch_Order_Parameters_b.(branch_parameter);
    
            % The analysed data
            data_labels         = {'mean rel. diff.', 'STD rel. diff.', 'total mean rel. diff.', 'total STD rel. diff.'};
            data_fields_cell    = {'MC.Mean.mean', 'MC.Mean.std', 'MC.Total.mean', 'MC.Total.std'};
            number_labels       = length(data_labels);
    
            relative_diff_cell  = cell(1, number_scanner_combinations);
    
            for c = 1 : number_scanner_combinations                
                % Data of each field for each branch order
                data_matrix_a = zeros(max_branch_order + 1, number_labels);
                data_matrix_b = zeros(max_branch_order + 1, number_labels);
    
                for b = 1 : max_branch_order
                    % This branch order's structures
                    Branch_Order_Structure_a = Branch_Parameter_Data_a(c, b);
                    Branch_Order_Structure_b = Branch_Parameter_Data_b(c, b);
    
                    % Data for each column
                    for f = 1 : number_labels
                        field               = data_fields_cell{f};
    
                        field_data_a        = Nested_Structure_Access(Branch_Order_Structure_a, field);
                        data_matrix_a(b, f) = field_data_a;
    
                        field_data_b        = Nested_Structure_Access(Branch_Order_Structure_b, field);
                        data_matrix_b(b, f) = field_data_b;
                    end
                end
    
                % The relative difference
                rel_diff_matrix         = 2*(data_matrix_a - data_matrix_b) ./ (data_matrix_a + data_matrix_b) * 100;
                relative_diff_cell{c}   = rel_diff_matrix;
    
                %--% Table %--%
                if Print == true
                    fprintf('%s [%s] \n', branch_parameter, branch_parameter_unit);
                end
    
                scanner_combination = scanner_combinations_cell{c};
                table_file_name = sprintf('TreeQSM_Rel_Diff_%s_branch_%s.xls', num2str(scanner_combination), branch_parameter);
    
                Table_Formatter(rel_diff_matrix, num_decimal_digits, output_format, table_row_labels, data_labels, table_file_name, Print);
        
                movefile(table_file_name, save_folder_name);
            end
    
            %--% Plots %--%
            if Plot == true
                for i = 1 : number_labels
                    % The data
                    data_label              = data_labels{i};
                    data_rel_diff_matrix    = zeros(max_branch_order + 1, number_scanner_combinations);
    
                    for c = 1 : number_scanner_combinations                
                        rel_diff_matrix             = relative_diff_cell{c};
                        data_rel_diff_matrix(:, c)  = rel_diff_matrix(:, i);
                    end
    
                    for s = 1 : number_branch_sets
                        % Plotted branch orders
                        branch_orders = branch_order_cell{s};
    
                        % Figure
                        figure_number = 100*p + 10*i + s;
                        figure(figure_number)
                        % Set the size and white background color
                        set(gcf, 'Units', 'Normalized', 'Position', [0.0 0.05 0.95 0.8])
                        set(gcf, 'color', [1, 1, 1])  
                        
                        hold on
                        grid on
    
                        % Bar graph
                        Bar_Graph = bar(branch_orders, data_rel_diff_matrix(branch_orders + 1, :));
    
                        for c = 1 : number_scanner_combinations
                            Bar_Graph(c).DisplayName    = scanner_combination_labels{c};
                            Bar_Graph(c).FaceColor      = scanner_combination_cmap(c, :);
                        end
    
                        % Axes
                        max_branch_order = max(branch_orders);
                        min_branch_order = min(branch_orders);
                        xlim([min_branch_order - 0.5, max_branch_order + 0.5]);
    
                        if max_branch_order == 0
                            xlabel('stem');
                            xticklabels('');
                            xticks('');
                        else
                            xlabel('branch order [-]');
                        end
    
                        branch_parameter_label = regexprep(branch_parameter, '([A-Z])', ' ${lower($1)}');           % Inserts a space before each capital letter and makes them lower case
                        
                        if strcmp(branch_parameter_label(1), ' ')         % This may result in a space being placed at the start, which is removed
                            branch_parameter_label(1) = [];
                        end
        
                        ylabel(sprintf('%s %s [%%]', branch_parameter_label, data_label));
        
                        % Legend
                        legend('show', 'location', 'eastoutside');
        
                        set(gca, 'FontSize', 15);
                        set(gca, 'LineWidth', 2);
        
                        hold off
        
                        % Saving the figure
                        if max_branch_order == 0
                            figure_name = sprintf('TreeQSM_stem_%s_%s', branch_parameter, data_label);
                        else
                            figure_name = sprintf('TreeQSM_branch_%s_%s', branch_parameter, data_label);
                        end
        
                        if c == number_scanner_combinations
                            export_fig(figure_number, [figure_name, '.fig']);
                            export_fig(figure_number, [figure_name, '.png']);
            
                            movefile([figure_name, '*'], save_folder_name); 
                        end
                    end
                end
            end
        end

    %% Number branches %%
        %--% Tables %--%
        [MC_number_branches_mean_matrix_a, MC_number_branches_std_matrix_a] = deal(Number_Branches_a.MC.mean, Number_Branches_a.MC.std);
        MC_number_branches_mean_matrix_a    = [MC_number_branches_mean_matrix_a, zeros(number_scanner_combinations, max_branch_order - max_branch_order_a)];
        MC_number_branches_std_matrix_a     = [MC_number_branches_std_matrix_a, NaN(number_scanner_combinations, max_branch_order - max_branch_order_a)];

        [MC_number_branches_mean_matrix_b, MC_number_branches_std_matrix_b] = deal(Number_Branches_b.MC.mean, Number_Branches_b.MC.std);
        MC_number_branches_mean_matrix_b    = [MC_number_branches_mean_matrix_b, zeros(number_scanner_combinations, max_branch_order - max_branch_order_b)];
        MC_number_branches_std_matrix_b     = [MC_number_branches_std_matrix_b, NaN(number_scanner_combinations, max_branch_order - max_branch_order_b)];

        % Mean relative difference
        MC_number_branches_mean_rel_diff_matrix = 2*(MC_number_branches_mean_matrix_a - MC_number_branches_mean_matrix_b) ./ (MC_number_branches_mean_matrix_a + MC_number_branches_mean_matrix_b) * 100;

        branch_order_labels = string(0:max_branch_order);
        table_file_name     = 'Number_Branches_MC_Mean_Rel_Diff.xls';
        Table_Formatter(MC_number_branches_mean_rel_diff_matrix, num_decimal_digits, output_format, scanner_combination_labels, branch_order_labels, table_file_name, Print);

        % Standard deviation
        MC_number_branches_std_rel_diff_matrix = 2*(MC_number_branches_std_matrix_a - MC_number_branches_std_matrix_b) ./ (MC_number_branches_std_matrix_a + MC_number_branches_std_matrix_b) * 100;

        table_file_name     = 'Number_Branches_MC_STD_Rel_Diff.xls';
        Table_Formatter(MC_number_branches_std_rel_diff_matrix, num_decimal_digits, output_format, scanner_combination_labels, branch_order_labels, table_file_name, Print);

        movefile('Number_Branches*Rel_Diff.xls', save_folder_name);

    %% Tree parameters %%
        % Various properties of the tree
        tree_parameter_fields   = {'TotalVolume', 'TrunkVolume', 'BranchVolume', 'TreeHeight', 'TrunkLength', 'BranchLength', 'TotalLength', 'NumberBranches', 'TrunkArea', 'BranchArea', 'TotalArea', 'DBHqsm', 'CrownDiamAve', 'CrownAreaConv', 'CrownBaseHeight', 'CrownLength'};
        tree_parameter_units    = {'m^3', 'm^3', 'm^3', 'm', 'm', 'm', 'm', '-', 'm^2', 'm^2', 'm^2', 'm', 'm', 'm^2', 'm', 'm'};
        number_tree_parameters  = length(tree_parameter_units);
    
        %--% Data %--%
        MC_tree_parameter_cell_a = cell(number_scanner_combinations, number_tree_parameters);
        MC_tree_parameter_cell_b = cell(number_scanner_combinations, number_tree_parameters);
    
        for t = 1 : number_tree_parameters
            % The parameter
            tree_parameter      = tree_parameter_fields{t};
            Tree_Parameter_fun  = @(QSM) QSM.treedata.(tree_parameter);
    
            for c = 1 : number_scanner_combinations    
                % Monte Carlo QSM data
                MC_QSM_cell_c_a                 = MC_QSM_cell_a{c};
                MC_tree_parameter_list_a        = cellfun(Tree_Parameter_fun, MC_QSM_cell_c_a);
                MC_tree_parameter_cell_a{c, t}  = MC_tree_parameter_list_a;
    
                MC_QSM_cell_c_b                 = MC_QSM_cell_b{c};
                MC_tree_parameter_list_b        = cellfun(Tree_Parameter_fun, MC_QSM_cell_c_b);
                MC_tree_parameter_cell_b{c, t}  = MC_tree_parameter_list_b;
            end
        end
    
        % Mean values and standard deviations
        MC_tree_parameter_mean_matrix_a         = cellfun(@mean, MC_tree_parameter_cell_a);
        MC_tree_parameter_mean_matrix_b         = cellfun(@mean, MC_tree_parameter_cell_b);
        MC_tree_parameter_mean_rel_diff_matrix  = 2*(MC_tree_parameter_mean_matrix_a - MC_tree_parameter_mean_matrix_b) ./ (MC_tree_parameter_mean_matrix_a + MC_tree_parameter_mean_matrix_b) * 100;
    
        MC_tree_parameter_std_matrix_a          = cellfun(@std, MC_tree_parameter_cell_a);
        MC_tree_parameter_std_matrix_b          = cellfun(@std, MC_tree_parameter_cell_b);
        MC_tree_parameter_std_rel_diff_matrix   = 2*(MC_tree_parameter_std_matrix_a - MC_tree_parameter_std_matrix_b) ./ (MC_tree_parameter_std_matrix_a + MC_tree_parameter_std_matrix_b) * 100;
    
        %--% Plots %--%
        if Plot == true
            for t = 1 : number_tree_parameters
                % The parameter and its data
                tree_parameter = tree_parameter_fields{t};
    
                % Figure                
                figure_number = 1000 + t;
                figure(figure_number)
                
                % Title
                tree_parameter_label = regexprep(tree_parameter, '([A-Z])', ' ${lower($1)}');           % Inserts a space before each capital letter and makes them lower case
                        
                if strcmp(tree_parameter_label(1), ' ')         % This may result in a space being placed at the start, which is removed
                    tree_parameter_label(1) = [];
                end

                sgtitle(tree_parameter_label);

                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])   
    
                % Mean values
                subplot(2, 1, 1)
                hold on
                grid on
    
                for c = 1 : number_scanner_combinations
                    scanner_comb_colour = scanner_combination_cmap(c, :);
                    bar(c, MC_tree_parameter_mean_rel_diff_matrix(c, t), 'FaceColor', scanner_comb_colour);
                end
    
                % Axes
                xticks(1 : number_scanner_combinations);
                xticklabels(scanner_combination_labels);
                ylabel('mean rel. diff. [%%]');
    
                % Text
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
    
                hold off
    
                % Standard deviations
                subplot(2, 1, 2)
                hold on
                grid on
    
                for c = 1 : number_scanner_combinations
                    scanner_comb_colour = scanner_combination_cmap(c, :);
                    bar(c, MC_tree_parameter_std_rel_diff_matrix(c, t), 'FaceColor', scanner_comb_colour);
                end
    
                % Axes
                xticks(1 : number_scanner_combinations);
                xticklabels(scanner_combination_labels);
                ylabel('std rel. diff. [%%]');
    
                % Text
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
    
                hold off
    
                % The figure is saved as .png and .fig and moved to the folder
                figure_name = sprintf('TreeQSM_tree_%s', tree_parameter);
                export_fig(figure_number, [figure_name, '.fig']);
                export_fig(figure_number, [figure_name, '.png']);
            
                movefile([figure_name, '*'], save_folder_name); 
            end
        end
    
        %--% Tables %--%
        % Mean absolute values
        Average_mean_rel_diff_list  = mean(abs(MC_tree_parameter_mean_rel_diff_matrix), 2);
        Average_std_rel_diff_list   = mean(abs(MC_tree_parameter_std_rel_diff_matrix), 2);
        
        table_file_name     = 'Tree_Parameter_Mean_Absolute.xls';
        average_table_data  = [Average_mean_rel_diff_list, Average_std_rel_diff_list];
        column_names        = {'MC mean', 'MC std'};

        Table_Formatter(average_table_data, num_decimal_digits, output_format, scanner_combination_labels, column_names, table_file_name, Print);
        movefile(table_file_name, save_folder_name);

        % Scanner combination data
        for c = 1 : number_scanner_combinations
            % Data 
            scanner_combination                     = scanner_combinations_cell{c};
            MC_tree_parameter_rel_diff_mean_list    = MC_tree_parameter_mean_rel_diff_matrix(c, :)';
            MC_tree_parameter_rel_diff_std_list     = MC_tree_parameter_std_rel_diff_matrix(c, :)';
            average_mean_rel_diff                   = Average_mean_rel_diff_list(c);
            average_std_rel_diff                    = Average_std_rel_diff_list(c);
    
            table_file_name     = sprintf('Tree_Parameter_%s_Results.xls', num2str(scanner_combination));
            tree_table_data     = [MC_tree_parameter_rel_diff_mean_list, MC_tree_parameter_rel_diff_std_list; average_mean_rel_diff, average_std_rel_diff];
            
            Table_Formatter(tree_table_data, num_decimal_digits, output_format, [tree_parameter_fields, 'mean abs.'], column_names, table_file_name, Print);
    
            movefile(table_file_name, save_folder_name);
        end

    %% Height-based results %%
        % Tree height of each QSM
        Tree_Height_fun     = @(QSM) QSM.treedata.TreeHeight;
        tree_height_cell_a  = cell(1, number_scanner_combinations);
        tree_height_cell_b  = cell(1, number_scanner_combinations);
    
        for c = 1 : number_scanner_combinations
            % This scanner combination's heights
            MC_QSM_cell_c_a         = MC_QSM_cell_a{c};
            tree_height_list_a      = cellfun(Tree_Height_fun, MC_QSM_cell_c_a);
            tree_height_cell_a{c}   = tree_height_list_a;
    
            MC_QSM_cell_c_b         = MC_QSM_cell_b{c};
            tree_height_list_b      = cellfun(Tree_Height_fun, MC_QSM_cell_c_b);
            tree_height_cell_b{c}   = tree_height_list_b;
        end
    
        tree_height_cell = [tree_height_cell_a; tree_height_cell_b];
        tree_height_list = horzcat(tree_height_cell{:});
        
        % Bins
        max_height          = max(tree_height_list);
        bin_edge_heights    = 0 : height_step : max_height;
    
        %--% The data for each scanner combination %--%
        QSM_Height_fun  = @(QSM) QSM_Height_Volume(QSM, bin_edge_heights);
        QSM_Height_Data = struct('branch_volume', [], 'branch_number_cylinders', [], 'stem_volume', [], 'stem_number_cylinders', []);
    
        for r = 1 : 2
            for c = 1 : number_scanner_combinations
                % This scanner combination's data
                if r == 1
                    MC_QSM_cell_c = MC_QSM_cell_a{c};
                elseif r == 2
                    MC_QSM_cell_c = MC_QSM_cell_b{c};
                end
        
                % Volume and number cylinders per height bin
                [branch_height_volume_cell, branch_height_num_cyl_cell, stem_height_volume_cell, stem_height_num_cyl_cell] = cellfun(QSM_Height_fun, MC_QSM_cell_c, 'UniformOutput', false);
                
                % Branches
                branch_height_volume_matrix     = vertcat(branch_height_volume_cell{:});
                branch_height_volume_mean_list  = mean(branch_height_volume_matrix, 1);
                branch_height_volume_std_list   = std(branch_height_volume_matrix, [], 1);
        
                numerical_margin                                = 1e-6;
                branch_height_volume_rel_uncert_list            = branch_height_volume_std_list ./ branch_height_volume_mean_list * 100;
                zero_bool                                       = branch_height_volume_mean_list < numerical_margin;                               
                branch_height_volume_rel_uncert_list(zero_bool) = 0;
        
                branch_height_num_cyl_matrix       = vertcat(branch_height_num_cyl_cell{:});
                branch_height_num_cyl_mean_list    = mean(branch_height_num_cyl_matrix, 1);
                branch_height_num_cyl_std_list     = std(branch_height_num_cyl_matrix, [], 1);
        
                branch_height_num_cyl_rel_uncert_list               = branch_height_num_cyl_std_list ./ branch_height_num_cyl_mean_list * 100;
                zero_bool                                           = branch_height_num_cyl_mean_list < numerical_margin;                               
                branch_height_num_cyl_rel_uncert_list(zero_bool)    = 0;
        
                % Stem
                stem_height_volume_matrix       = vertcat(stem_height_volume_cell{:});
                stem_height_volume_mean_list    = mean(stem_height_volume_matrix, 1);
                stem_height_volume_std_list     = std(stem_height_volume_matrix, [], 1);
        
                stem_height_volume_rel_uncert_list              = stem_height_volume_std_list ./ stem_height_volume_mean_list * 100;
                zero_bool                                       = stem_height_volume_mean_list < numerical_margin;                               
                stem_height_volume_rel_uncert_list(zero_bool)   = 0;
        
                stem_height_num_cyl_matrix 	    = vertcat(stem_height_num_cyl_cell{:});
                stem_height_num_cyl_mean_list   = mean(stem_height_num_cyl_matrix, 1);
                stem_height_num_cyl_std_list    = std(stem_height_num_cyl_matrix, [], 1);
        
                stem_height_num_cyl_rel_uncert_list             = stem_height_num_cyl_std_list ./ stem_height_num_cyl_mean_list * 100;
                zero_bool                                       = stem_height_num_cyl_mean_list < numerical_margin;                               
                stem_height_num_cyl_rel_uncert_list(zero_bool)  = 0;
        
        
                QSM_Height_Data(r, c) = struct('branch_volume', struct('mean', branch_height_volume_mean_list, 'std', branch_height_volume_std_list, 'relative_uncertainty', branch_height_volume_rel_uncert_list), 'branch_number_cylinders', struct('mean', branch_height_num_cyl_mean_list, 'std', branch_height_num_cyl_std_list, 'relative_uncertainty', branch_height_num_cyl_rel_uncert_list), ... 
                                               'stem_volume', struct('mean', stem_height_volume_mean_list, 'std', stem_height_volume_std_list, 'relative_uncertainty', stem_height_volume_rel_uncert_list), 'stem_number_cylinders', struct('mean', stem_height_num_cyl_mean_list, 'std', stem_height_num_cyl_std_list, 'relative_uncertainty', stem_height_num_cyl_rel_uncert_list));
            end
        end
    
        %--% Plots %--%
        if Plot == true
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
                data_fields         = fieldnames(QSM_Height_Data(1).(height_metric));
                number_data_fields  = length(data_fields);
    
                for d = 1 : number_data_fields
                    % Plotting the data
                    data_field = data_fields{d};
    
                    %--% Absolute values %--%
                    figure_number = 2000 + 100*m + d*10;
                    figure(figure_number)
    
                    % Set the size and white background color
                    set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.9])
                    set(gcf, 'color', [1, 1, 1])   
    
                    Tiled_Chart_Layout_Main = tiledlayout(1, 1);
                    Tiled_Chart_Layout      = gobjects(1, 1);
                    Axes                    = gobjects([1, number_scanner_combinations]);
    
                    Tiled_Chart_Layout(1) = tiledlayout(Tiled_Chart_Layout_Main, 1, number_scanner_combinations);
                    Tiled_Chart_Layout(1).Layout.Tile = 1;
    
                    % Each scanner combination's and folder's data
                    for c = 1 : number_scanner_combinations
                        scanner_combination_label = scanner_combination_labels{c};
    
                        Axes(1, c) = nexttile(Tiled_Chart_Layout(1));
                        hold on
                        grid on
    
                        title(scanner_combination_label)
    
                        for r = 1 : 2
                            results_label = QSM_labels_cell{r};
    
                            % The data
                            height_data_list = QSM_Height_Data(r, c).(height_metric).(data_field);
        
                            % Plot
                            plot(height_data_list, bin_avg_heights, 'LineWidth', 2, 'DisplayName', results_label);
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
            
                        % Legend and text
                        if c == 1
                            ylabel('H [m]');
                            legend('show', 'location', 'northeast');
                        end
        
                        hold off    
                    end
    
                    % Shared x-axis label
                    x_label = sprintf('%s %s [%s]', metric_label, field_label, data_unit);
                    xlabel(Tiled_Chart_Layout(1), x_label, 'FontSize', 15);
    
                    % Formatting
                    set(Axes, 'FontSize', 15);
                    set(Axes, 'LineWidth', 2);
    
                    % Saving the figure
                    figure_name = sprintf('TreeQSM_Height_%s_%s', height_metric, data_field);
                    export_fig(figure_number, [figure_name, '.fig']);
                    export_fig(figure_number, [figure_name, '.png']);
    
                    movefile([figure_name, '*'], save_folder_name); 
    
                    %--% Relative difference %--%
                    figure_number = 3000 + 100*m + d*10;
                    figure(figure_number)
    
                    % Set the size and white background color
                    set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.9])
                    set(gcf, 'color', [1, 1, 1])   
    
                    Tiled_Chart_Layout_Main = tiledlayout(1, 1);
                    Tiled_Chart_Layout      = gobjects(1, 1);
                    Axes                    = gobjects([1, number_scanner_combinations]);
    
                    Tiled_Chart_Layout(1) = tiledlayout(Tiled_Chart_Layout_Main, 1, number_scanner_combinations);
                    Tiled_Chart_Layout(1).Layout.Tile = 1;
    
                    % Each scanner combination's and folder's data
                    for c = 1 : number_scanner_combinations
                        % Scanner combination
                        scanner_combination_label   = scanner_combination_labels{c};
                        scanner_comb_colour         = scanner_combination_cmap(c, :);
    
                        % Relative difference
                        height_data_list_a = QSM_Height_Data(1, c).(height_metric).(data_field);
                        height_data_list_b = QSM_Height_Data(2, c).(height_metric).(data_field);
    
                        height_data_rel_diff_list = 2*(height_data_list_a - height_data_list_b) ./ (height_data_list_a + height_data_list_b) * 100;
        
                        % Subplot
                        Axes(1, c) = nexttile(Tiled_Chart_Layout(1));
                        hold on
                        grid on
    
                        title(scanner_combination_label)
        
                        % Plot
                        plot(height_data_rel_diff_list, bin_avg_heights, 'LineWidth', 2, 'Color', scanner_comb_colour);
    
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
            
                        if c == 1
                            ylabel('H [m]');
                        end
        
                        hold off    
                    end
    
                    % Shared x-axis label
                    x_label = sprintf('%s %s \n rel. diff. [%%]', metric_label, field_label);
                    xlabel(Tiled_Chart_Layout(1), x_label, 'FontSize', 15);
    
                    % Formatting
                    set(Axes, 'FontSize', 15);
                    set(Axes, 'LineWidth', 2);
    
                    % Saving the figure
                    figure_name = sprintf('TreeQSM_Height_%s_%s_Rel_Diff', height_metric, data_field);
                    export_fig(figure_number, [figure_name, '.fig']);
                    export_fig(figure_number, [figure_name, '.png']);
    
                    movefile([figure_name, '*'], save_folder_name); 
                end
            end
        end
end