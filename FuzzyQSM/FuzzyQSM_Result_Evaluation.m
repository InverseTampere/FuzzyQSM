% This script shows the initial discrete QSM and fuzzy QSM along with the point cloud, 
% as well as how it compares to the original model

function Results_Tables = FuzzyQSM_Result_Evaluation(tree_ID, Original_QSM, Discrete_QSM, Fuzzy_QSM, Smoothed_Model_Metrics, Smoothed_Model_Circle_Geometry, Cyl_Point_Cloud_Distributions_cell, Data_Parameters, Fitting_Parameters, Scanning_Parameters, Scanner_Parameters, Parallel_Pool, Output_Decisions)

    %% Structure inputs %%
        % Smoothed model
        Smoothed_Model_Tree_Metrics         = Smoothed_Model_Metrics.treedata;
        Smoothed_Model_Branch_Order_Metrics = Smoothed_Model_Metrics.BranchOrder;

        % Discrete QSM
        Discrete_Cylinder_Geometry          = Discrete_QSM.cylinder;

        % Fuzzy QSM
        Fuzzy_Cylinder_Geometry             = Fuzzy_QSM.cylinder;

        % Outputs
        Plot                                = Output_Decisions.Plot;
        Print                               = Output_Decisions.Print;
        Compute_Radius_Errors               = Output_Decisions.Compute_Radius_Errors;

    %% Manual inputs %%
        % Point cloud settings
        coarsening_factor           = 1e1;          % [-] Factor by which the point cloud is coarsened. Recommended to be fairly high

        % Table settings
        num_decimal_digits          = 8;            % [-]
        output_format               = 'Float';      % [Float, Integer, Exponential]

        % Height binning
        max_height                  = [];         % [m] Optional. Maximum height used for binning
        height_step                 = 0.2;          % [m] Size of each bin
        confidence_interval         = 90;           % [%] Shows the standard deviation for the geometric uncertainty

        % Box plots
        Separate_Stem               = true;        % [true, false] Separate the stem from the branches for the plots

    %% Preliminary %%
        % A folder is created with the current time as its name to save the results in
        current_time        = datetime('now', 'format', 'yyyyMMdd_HHmmss');                 % Date followed by clock time       
        save_folder_name    = sprintf('%s_%s', tree_ID, string(current_time));              % Converted to string array that acts as the name of the folder

        mkdir(save_folder_name);

        % Inputs of the cylinder fitting algorithm are saved
        save('FuzzyQSM_Inputs', 'tree_ID', 'Data_Parameters', 'Fitting_Parameters', 'Scanning_Parameters', 'Scanner_Parameters');
        movefile('FuzzyQSM_Inputs.mat', save_folder_name);

    %% QSM metrics %%
        % Extended branch and tree metrics for the QSMs
        [~, Discrete_Branch_Order_Metrics]  = QSM_Branch_Metrics(Discrete_Cylinder_Geometry, Discrete_QSM);
        Discrete_Tree_Metrics               = QSM_Tree_Metrics(Discrete_Cylinder_Geometry, Discrete_Branch_Order_Metrics);

        Fuzzy_QSM.branch.order              = Discrete_QSM.branch.order;
        [~, Fuzzy_Branch_Order_Metrics]     = QSM_Branch_Metrics(Fuzzy_Cylinder_Geometry, Fuzzy_QSM);
        Fuzzy_Tree_Metrics                  = QSM_Tree_Metrics(Fuzzy_Cylinder_Geometry, Fuzzy_Branch_Order_Metrics);

        % Maximum branch order
        max_discrete_branch_order   = Discrete_Branch_Order_Metrics.max_branch_order;
        max_fuzzy_branch_order      = Fuzzy_Branch_Order_Metrics.max_branch_order;
        max_obj_branch_order        = Smoothed_Model_Branch_Order_Metrics.max_branch_order;
        max_branch_order            = max([max_discrete_branch_order, max_fuzzy_branch_order, max_obj_branch_order]);

    %% QSMs %%
        % The analysed QSMs
        QSM_cell    = {Original_QSM, Discrete_QSM, Fuzzy_QSM};
        QSM_labels  = {'Original', 'Discrete', 'Fuzzy'};
        number_QSMs = length(QSM_labels);

        % The point cloud
        Distribution_Mu_fun     = @(Cyl_Point_Cloud_Distributions) vertcat(Cyl_Point_Cloud_Distributions.distribution_mu_cell{:});
        empty_bool              = cellfun(@isempty, Cyl_Point_Cloud_Distributions_cell);
        distribution_mu_cell    = cellfun(Distribution_Mu_fun, Cyl_Point_Cloud_Distributions_cell(~empty_bool), 'UniformOutput', false);
        distribution_mu_matrix  = vertcat(distribution_mu_cell{:});

        if Plot == true
            % QSM plot inputs
            figure_number   = 1;
            alpha_value     = 1.0;
            number_facets   = 10;
           
            % Figure
            Fig         = figure(figure_number);
            figure_name = 'Fitted_QSMs';
  
            set(Fig, 'name', figure_name, 'NumberTitle', 'off');

            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    

            for q = 1 : number_QSMs
                % Plotted QSM
                QSM         = QSM_cell{q};
                QSM_label   = QSM_labels{q};

                subplot(1, number_QSMs, q)
                hold on
                grid on
                title(QSM_label);
        
                % QSM
                plot_cylinder_model(QSM.cylinder, "order", figure_number, number_facets, alpha_value);

                % The point cloud
                scatter3(distribution_mu_matrix(1:coarsening_factor:end, 1), distribution_mu_matrix(1:coarsening_factor:end, 2), distribution_mu_matrix(1:coarsening_factor:end, 3), 1, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');

                % Axes
                xlabel('x [m]');
                ylabel('y [m]');
                zlabel('z [m]');

                view(45, 45);

                % Formatting
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off
            end

            export_fig(figure_number, [figure_name, '.fig']);
            export_fig(figure_number, [figure_name, '.png']);

            movefile([figure_name, '*'], save_folder_name); 
        end

    %% Tree metric evaluation %%
        % Tree metrics
        tree_metrics        = fieldnames(Smoothed_Model_Tree_Metrics);
        number_tree_metrics = length(tree_metrics);

        [discrete_tree_metrics_list, fuzzy_tree_metrics_list, qsm_obj_tree_metrics_list] = deal(zeros(1, number_tree_metrics));

        for t = 1 : number_tree_metrics
            tree_metric                     = tree_metrics{t};
            discrete_tree_metrics_list(t)   = Discrete_Tree_Metrics.(tree_metric);
            fuzzy_tree_metrics_list(t)      = Fuzzy_Tree_Metrics.(tree_metric);
            qsm_obj_tree_metrics_list(t)    = Smoothed_Model_Tree_Metrics.(tree_metric);
        end

        table_row_labels    = {'Discrete', 'Fuzzy', 'Original'};
        tree_metrics_matrix = [discrete_tree_metrics_list; fuzzy_tree_metrics_list; qsm_obj_tree_metrics_list];
        table_file_name     = 'Tree_Metrics.xls';
        Tree_Metrics_Table  = Table_Formatter(tree_metrics_matrix, num_decimal_digits, output_format, table_row_labels, tree_metrics, table_file_name, Print);
        
        Results_Tables.Tree_Metrics.Values = Tree_Metrics_Table;
        movefile(table_file_name, save_folder_name);

        % Relative errors
        discrete_tree_metrics_rel_error_list    = (discrete_tree_metrics_list - qsm_obj_tree_metrics_list) ./ qsm_obj_tree_metrics_list * 100;
        fuzzy_tree_metrics_rel_error_list       = (fuzzy_tree_metrics_list - qsm_obj_tree_metrics_list) ./ qsm_obj_tree_metrics_list * 100;

        table_row_labels                = {'Discrete', 'Fuzzy'};
        tree_metrics_rel_error_matrix   = [discrete_tree_metrics_rel_error_list; fuzzy_tree_metrics_rel_error_list];
        table_file_name                 = 'Tree_Metrics_Relative_Errors.xls';

        Tree_Metrics_Rel_Error_Table    = Table_Formatter(tree_metrics_rel_error_matrix, num_decimal_digits, output_format, table_row_labels, tree_metrics, table_file_name, Print);

        Results_Tables.Tree_Metrics.Rel_Error = Tree_Metrics_Rel_Error_Table;
        movefile(table_file_name, save_folder_name);

    %% Branch-order metric evaluation %%
        % Branch orders
        branch_order_labels = string(0 : max_branch_order);

        % Evaluated branch order metrics
        branch_order_metrics        = fieldnames(Discrete_Branch_Order_Metrics);
        number_branch_order_metrics = length(branch_order_metrics);

        for m = 1 : number_branch_order_metrics
            % Assessed metric
            branch_order_metric = branch_order_metrics{m};

            if strcmp(branch_order_metric, 'max_branch_order')          % Skip the maximum branch order field
                continue
            end

            % Metric values
            discrete_metric_list    = zeros(max_branch_order + 1, 1);
            fuzzy_metric_list       = zeros(max_branch_order + 1, 1);
            qsm_obj_metric_list     = zeros(max_branch_order + 1, 1);

            for b = 1 : max_branch_order + 1
                if b <= Discrete_Branch_Order_Metrics.max_branch_order + 1
                    discrete_metric_list(b) = Discrete_Branch_Order_Metrics.(branch_order_metric)(b);
                end
                if b <= Fuzzy_Branch_Order_Metrics.max_branch_order + 1
                    fuzzy_metric_list(b)    = Fuzzy_Branch_Order_Metrics.(branch_order_metric)(b);
                end
                if b <= Smoothed_Model_Branch_Order_Metrics.max_branch_order + 1
                    qsm_obj_metric_list(b)  = Smoothed_Model_Branch_Order_Metrics.(branch_order_metric)(b);
                end
            end

            % Relative errors
            discrete_rel_error_list = (discrete_metric_list - qsm_obj_metric_list) ./ qsm_obj_metric_list * 100;
            fuzzy_rel_error_list    = (fuzzy_metric_list - qsm_obj_metric_list) ./ qsm_obj_metric_list * 100;

            % Table
            table_column_labels = {'Original', 'Discrete', 'Fuzzy', 'Discr. rel. error', 'Fuzzy rel. error'};
            table_matrix        = [qsm_obj_metric_list, discrete_metric_list, fuzzy_metric_list, discrete_rel_error_list, fuzzy_rel_error_list];
            table_file_name     = sprintf('Branch_Order_%s.xls', branch_order_metric);

            Metrics_Table = Table_Formatter(table_matrix, num_decimal_digits, output_format, branch_order_labels, table_column_labels, table_file_name, Print);

            Results_Tables.Branch_Order_Metrics.(branch_order_metric) = Metrics_Table;
            movefile(table_file_name, save_folder_name);
        end

    %% Fitting height histograms %%        
        % Bins
        if isempty(max_height)
            % If not provided, the maximum height is taken from the QSMs
            discrete_tree_height    = Discrete_QSM.treedata.TreeHeight;
            fuzzy_tree_height       = Fuzzy_QSM.treedata.TreeHeight;
            original_tree_height    = Original_QSM.treedata.TreeHeight;

            max_height = max([discrete_tree_height, fuzzy_tree_height, original_tree_height]);
        end

        bin_edge_heights = 0 : height_step : max_height;

        % Data for each QSM 
        QSM_Height_Data = struct();

        for q = 1 : number_QSMs
            % This QSM's data in each bin
            QSM         = QSM_cell{q};
            QSM_label   = QSM_labels{q};

            [branch_total_volume_list, branch_number_cylinders_list, stem_total_volume_list, stem_number_cylinders_list] = QSM_Height_Volume(QSM, bin_edge_heights);

            QSM_Height_Data.(QSM_label) = struct('branch_total_volume', branch_total_volume_list, 'branch_number_cylinders', branch_number_cylinders_list, 'stem_total_volume', stem_total_volume_list, 'stem_number_cylinders', stem_number_cylinders_list);
        end

        % Saved in a structure
        save('QSM_Height_Data.mat', 'QSM_Height_Data', 'bin_edge_heights');
        movefile('QSM_Height_Data.mat', save_folder_name);
    
        %--% Plots %--%
        if Plot == true
            % Average height for each bin
            bin_avg_heights = (bin_edge_heights(2 : end) + bin_edge_heights(1 : end - 1)) / 2;

            % Colours for each QSM
            QSM_colour_map = cbrewer('qual', 'Set1', max(3, number_QSMs));
            QSM_colour_map = max(0, QSM_colour_map);
            QSM_colour_map = min(1, QSM_colour_map);

            % Plotted data
            height_metrics_cell = fieldnames(QSM_Height_Data.(QSM_label));
            number_metrics      = length(height_metrics_cell); 

            for m = 1 : number_metrics
                % The plotted height metric
                height_metric   = height_metrics_cell{m};
                metric_label    = strrep(height_metric, '_', ' ');     % Replace underscores with spaces

                if contains(metric_label, 'volume')
                    metric_unit = 'm^3';
                elseif contains(metric_label, 'number')
                    metric_unit = '-';
                else
                    warning('No unit specified for %s', metric_label);
                end

                % This metric's figure
                figure_number   = figure_number + 1;
                Fig             = figure(figure_number);
                figure_name     = sprintf('Height_%s', height_metric);

                set(Fig, 'name', figure_name, 'NumberTitle', 'off');

                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.5 0.9])
                set(gcf, 'color', [1, 1, 1])   

                hold on
                grid on

                % This metric's data for each QSM
                for q = 1 : number_QSMs
                    % QSM data
                    QSM_colour  = QSM_colour_map(q, :);
                    QSM_label   = QSM_labels{q};
                    QSM_data    = QSM_Height_Data.(QSM_label).(height_metric);

                    % Plot
                    plot(QSM_data, bin_avg_heights, 'LineWidth', 2, 'Color', QSM_colour, 'DisplayName', QSM_label);
                end

                x_axis_label = sprintf('%s [%s]', metric_label, metric_unit);
                xlabel(x_axis_label);
                set(gca, 'XTickLabelRotation', 0);

                ylim([0, max_height])
                ylabel('H [m]');

                % Legend and text
                legend('show', 'location', 'northoutside', 'orientation', 'horizontal');

                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off    

                % Saving the figure
                export_fig(figure_number, sprintf('%s.fig', figure_name));
                export_fig(figure_number, sprintf('%s.png', figure_name));

                movefile([figure_name, '*'], save_folder_name); 
            end
        end

    %% Uncertainty height histograms %%
        % Geometric uncertainty 
        Sigma_Radial_fun    = @(Cyl_Point_Cloud_Distributions) vertcat(Cyl_Point_Cloud_Distributions.sigma_radial_cell{:});
        sigma_radial_cell   = cellfun(Sigma_Radial_fun, Cyl_Point_Cloud_Distributions_cell(~empty_bool), 'UniformOutput', false);
        sigma_radial_list   = vertcat(sigma_radial_cell{:});

        Sigma_Prop_fun      = @(Cyl_Point_Cloud_Distributions) vertcat(Cyl_Point_Cloud_Distributions.sigma_prop_cell{:});
        sigma_prop_cell     = cellfun(Sigma_Prop_fun, Cyl_Point_Cloud_Distributions_cell(~empty_bool), 'UniformOutput', false);
        sigma_prop_list     = vertcat(sigma_prop_cell{:});

        geom_avg_sigma_list = sqrt(sigma_prop_list .* sigma_radial_list);

        % Binning in terms of height
        mu_height_list = distribution_mu_matrix(:, 3);

        number_bins                                     = length(bin_edge_heights) - 1;
        [geom_avg_sigma_mean_list, ~, bin_bool_cell]    = Data_Binning(geom_avg_sigma_list, mu_height_list, bin_edge_heights, number_bins);
        
        percentile_lower   = 50 - confidence_interval/2;
        percentile_upper   = 50 + confidence_interval/2;

        geom_avg_sigma_CI_matrix = zeros(number_bins, 2);
        
        for b = 1 : number_bins
            bin_bool                    = bin_bool_cell{b};
            geom_avg_sigma_bin_list     = geom_avg_sigma_list(bin_bool);

            sigma_confidence_interval       = prctile(geom_avg_sigma_bin_list, [percentile_lower, percentile_upper]);
            geom_avg_sigma_CI_matrix(b, :)  = sigma_confidence_interval;
        end

        % Save in a structure
        save('Height_Uncertainty.mat', 'sigma_radial_list', 'sigma_prop_list', 'geom_avg_sigma_list', 'mu_height_list');
        movefile('Height_Uncertainty.mat', save_folder_name);

        % Plot
        if Plot == true
            figure_number   = figure_number + 1;
            Fig             = figure(figure_number);
            figure_name     = 'Geom_Avg_Sigma';
        
            set(Fig, 'name', figure_name, 'NumberTitle', 'off');
        
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.6 0.9])
            set(gcf, 'color', [1, 1, 1])   

            hold on
            grid on

            % Mean value
            plot(geom_avg_sigma_mean_list * 1e3, bin_avg_heights, 'LineWidth', 2, 'color', 'b', 'DisplayName', sprintf('%s', '\mu'));       % Conversion to mm

            % Standard deviation
            [height_sigma_lower_list, height_sigma_upper_list] = Column_Deal(1e3 * geom_avg_sigma_CI_matrix);                               % Conversion to mm
            NaN_bool = isnan(height_sigma_lower_list);
            CI_label = sprintf('%i%% CI', confidence_interval);
            patch('XData', [height_sigma_lower_list(~NaN_bool); flipud(height_sigma_upper_list(~NaN_bool))], 'YData', [bin_avg_heights(~NaN_bool), fliplr(bin_avg_heights(~NaN_bool))]', 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.5, 'DisplayName', CI_label);

            % Axes
            xlabel('$\bar{\sigma}$ [mm]', 'interpreter', 'latex');
            xlim([0, max(height_sigma_upper_list)]);                   
            
            ylabel('H [m]');
            ylim([0, max_height]);
        
            % Legend
            legend('show', 'location', 'northoutside', 'orientation', 'horizontal');
        
            % Formatting
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
        
            hold off    
        
            % Saving the figure
            figure_name = 'Height_Uncertainty';
            export_fig(figure_number, sprintf('%s.fig', figure_name));
            export_fig(figure_number, sprintf('%s.png', figure_name));

            movefile([figure_name, '*'], save_folder_name); 
        end
   
    %% Radius error %%
        % Difference in radius between the object circles and nearest QSM cylinders
        if Compute_Radius_Errors == true
            Radius_Error_Diagnostics = false;
            [Discr_Radius_Errors, ~] = QSM_Radius_Error(Discrete_QSM, Smoothed_Model_Circle_Geometry, Parallel_Pool, Radius_Error_Diagnostics);
            [Fuzzy_Radius_Errors, ~] = QSM_Radius_Error(Fuzzy_QSM, Smoothed_Model_Circle_Geometry, Parallel_Pool, Radius_Error_Diagnostics);
    
            % Tables
            Radius_Errors_cell  = {Discr_Radius_Errors, Fuzzy_Radius_Errors};
            method_labels       = {'Discrete', 'Fuzzy'};
            number_methods      = length(method_labels);
    
            error_fields    = fieldnames(Discr_Radius_Errors.Total);
            number_errors   = length(error_fields);
    
            error_values    = {'mean', 'mean_abs'};
            number_values   = length(error_values);
    
            for e = 1 : number_errors
                % Each type of error (relative and regular) gets its own table
                error_field = error_fields{e};
    
                table_column_labels = cell(number_values, number_methods);
                table_errors_cell   = repmat({NaN(max_branch_order + 2, 1)}, [number_values, number_methods]);          % NaN values are specified for if a QSM does not have the maximum branch order
    
                for v = 1 : number_values
                    error_value = error_values{v};
    
                    for m = 1 : number_methods
                        method                      = method_labels{m};
                        column_label                = sprintf('%s %s radius %s', method, error_value, error_field);
                        table_column_labels{v, m}   = column_label;
    
                        Radius_Errors       = Radius_Errors_cell{m};
                        total_error         = Radius_Errors.Total.(error_field).(error_value);
                        
                        Branch_Order_cell   = {Radius_Errors.Branch_Order.(error_field)};
                        Error_Value_fun     = @(Errors) Errors.(error_value);
                        branch_order_errors = cellfun(Error_Value_fun, Branch_Order_cell);
    
                        error_list          = [total_error, branch_order_errors];
                        table_errors_cell{v, m}(1 : length(error_list)) = error_list;
                    end
                end
    
                % Table
                table_row_labels    = ["Total", branch_order_labels];
                table_column_labels = reshape(table_column_labels', [1, number_methods * number_values]);
                table_errors_cell   = reshape(table_errors_cell', [1, number_methods * number_values]);
                table_errors_matrix = horzcat(table_errors_cell{:});
                table_file_name     = sprintf('Radius_%s.xls', error_field);
                Radius_Error_Table  = Table_Formatter(table_errors_matrix, num_decimal_digits, output_format, table_row_labels, table_column_labels, table_file_name, Print);
                
                Results_Tables.Radius.(error_field) = Radius_Error_Table;
                movefile(table_file_name, save_folder_name);
            end
        end

        % Plot        
        if Plot == true && Compute_Radius_Errors == true
            %--% 3D error plot %--%
            % Colour maps used for the errors
            number_colours  = 1e3;
            error_cmap      = cbrewer('div', 'RdBu', number_colours);
            error_cmap      = max(error_cmap, 0);
            error_cmap      = min(error_cmap, 1);

            % Circle centres
            obj_circle_centre_matrix = Smoothed_Model_Circle_Geometry.centre;

            % Figure
            figure_number   = figure_number + 1;
            Fig             = figure(figure_number);
            figure_name     = 'Radius_Errors';
  
            set(Fig, 'name', figure_name, 'NumberTitle', 'off');

            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85]);
            set(gcf, 'color', [1, 1, 1]);

            % Subplot layout
            Tiled_Chart_Layout_Main = tiledlayout(number_methods, 1);
            Tiled_Chart_Layout      = gobjects(1, number_methods);
            Axes                    = gobjects([number_methods, number_errors]);

            % Subplots for each method and error type
            for m = 1 : number_methods
                % The method's row
                method_label = method_labels{m};

                Tiled_Chart_Layout(m) = tiledlayout(Tiled_Chart_Layout_Main, 1, number_errors);
                Tiled_Chart_Layout(m).Layout.Tile = m;
                
                for e = 1 : number_errors
                    % This error's subplot
                    Axes(m, e)  = nexttile(Tiled_Chart_Layout(m));

                    error_field = error_fields{e};
                    error_label = strrep(error_field, '_', ' ');

                    if contains(error_label, 'relative')
                        error_unit = '%';
                    else
                        error_unit = '-';
                    end

                    error_label = sprintf('%s [%s]', error_label, error_unit);

                    Error_fun   = @(Radius_Errors) Radius_Errors.Total.(error_field).data;
                    error_cell  = cellfun(Error_fun, Radius_Errors_cell, 'UniformOutput', false);
                    error_list  = error_cell{m};

                    NaN_bool             = isnan(error_list);
                    error_list(NaN_bool) = [];

                    error_list_total    = vertcat(error_cell{:});
                    error_limit         = max(abs(error_list_total));

                    error_list_n        = (error_list + error_limit) / (2*error_limit);
                    error_colour_ind    = round(error_list_n * (number_colours - 1)) + 1;
                    error_colours       = error_cmap(error_colour_ind, :);
        
                    hold on
                    grid on
        
                    sc_size = 1e1;
                    scatter3(obj_circle_centre_matrix(~NaN_bool, 1), obj_circle_centre_matrix(~NaN_bool, 2), obj_circle_centre_matrix(~NaN_bool, 3), sc_size, error_colours);
        
                    cb = colorbar;
                    shading interp
                    clim(error_limit * [-1, 1])
                    ylabel(cb, error_label);
                    cb.FontSize = 15;      
        
                    xlabel('x [m]');
                    ylabel('y [m]');
                    zlabel('z [m]');
        
                    axis equal
                    view(45, 45);                    
                end

                title(Tiled_Chart_Layout(m), method_label, 'FontSize', 15);
            end

            set(Axes, 'FontSize', 15);
            set(Axes, 'LineWidth', 2);

            colormap(error_cmap);

            % Save the figure
            export_fig(figure_number, [figure_name, '.fig']);
            export_fig(figure_number, [figure_name, '.png']);

            movefile([figure_name, '*'], save_folder_name); 

            %--% Branch-order box plots %--%
            % Box plot settings
            box_spacing     = 1/number_methods - 0.05;
            box_width       = box_spacing - 0.05;
            total_spacing   = number_methods * box_spacing;

            % Colours used for the methods
            method_colours = cbrewer('qual', 'Set1', max(3, number_methods));
            method_colours = max(method_colours, 0);
            method_colours = min(method_colours, 1);

            % Plotted branch orders
            if Separate_Stem == true
                branch_order_cell = {0, 1:max_branch_order};
            else
                branch_order_cell = {0:max_branch_order};
            end

            number_branch_sets = length(branch_order_cell);

            for s = 1 : number_branch_sets
                % Plotted branch orders
                branch_orders           = branch_order_cell{s};
                number_branch_orders    = length(branch_orders);

                % Figure
                figure_number   = figure_number + 1;
                Fig             = figure(figure_number);

                if max(branch_orders) == 0
                    figure_name = 'Radius_Errors_Stem';
                else
                    figure_name = 'Radius_Errors_Branch';
                end
      
                set(Fig, 'name', figure_name, 'NumberTitle', 'off');
    
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85]);
                set(gcf, 'color', [1, 1, 1]);
    
                hold on
                grid on
    
                for b = 1 : number_branch_orders
                    branch_order = branch_orders(b);
    
                    for m = 1 : number_methods
                        % The data
                        method_colour   = method_colours(m, :);
                        method_label    = method_labels{m};

                        try
                            error_list  = Radius_Errors_cell{m}.Branch_Order(branch_order + 1).error.data;
                        catch
                            error_list  = NaN;
                        end
    
                        % Location of the box plot
                        x_box = double(branch_order) - total_spacing/2 - box_spacing/2 + box_spacing * m; 
    
                        % Box plot
                        error_list = double(error_list);
                        box = boxchart(x_box * ones(length(error_list), 1), error_list, 'BoxFaceColor', method_colour, 'MarkerColor', method_colour, 'BoxWidth', box_width, 'Notch', 'off', 'DisplayName', method_label);
    
                        if b > 1
                            box.HandleVisibility = 'off';
                        end
                    end
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

                ylabel('Radius error [m]');

                % Legend
                legend('show', 'location', 'eastoutside');
    
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
    
                hold off
    
                % Saving the figure     
                export_fig(figure_number, sprintf('%s.fig', figure_name));
                export_fig(figure_number, sprintf('%s.png', figure_name));
    
                movefile([figure_name, '*'], save_folder_name); 
            end
        
            %--% Relationship with number of distributions %--%
            number_cylinders        = length(Cyl_Point_Cloud_Distributions_cell);
            cyl_number_points_list  = NaN(1, number_cylinders);
            empty_bool              = cellfun(@isempty, Cyl_Point_Cloud_Distributions_cell);
            
            Number_Distributions_fun            = @(Cyl_Point_Cloud_Distributions) Cyl_Point_Cloud_Distributions.number_distributions;
            cyl_number_points_list(~empty_bool) = cellfun(Number_Distributions_fun, Cyl_Point_Cloud_Distributions_cell(~empty_bool));

            % Figure
            figure_number   = figure_number + 1;
            Fig             = figure(figure_number);
            figure_name     = 'Radius_Error_Number_Distributions';
  
            set(Fig, 'name', figure_name, 'NumberTitle', 'off');

            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85]);
            set(gcf, 'color', [1, 1, 1]);

            % Subplots for each method and error type
            Tiled_Chart_Layout_Main = tiledlayout(number_methods, 1);
            Tiled_Chart_Layout      = gobjects(1, number_methods);
            Axes                    = gobjects([number_methods, number_errors]);

            for m = 1 : number_methods
                % The method's row
                method_label = method_labels{m};

                Tiled_Chart_Layout(m) = tiledlayout(Tiled_Chart_Layout_Main, 1, number_errors);
                Tiled_Chart_Layout(m).Layout.Tile = m;
                
                for e = 1 : number_errors
                    % This error's subplot
                    Axes(m, e)  = nexttile(Tiled_Chart_Layout(m));

                    error_field = error_fields{e};
                    error_label = strrep(error_field, '_', ' ');

                    if contains(error_label, 'relative')
                        error_unit = '%';
                    else
                        error_unit = '-';
                    end

                    error_label = sprintf('%s [%s]', error_label, error_unit);

                    Error_fun   = @(Radius_Errors) Radius_Errors.Cylinders.(error_field);
                    error_cell  = cellfun(Error_fun, Radius_Errors_cell, 'UniformOutput', false);
                    error_list  = error_cell{m};

                    error_list_total    = horzcat(error_cell{:});
                    error_max           = max(error_list_total);
                    error_min           = min(error_list_total);
        
                    hold on
                    grid on
        
                    scatter(cyl_number_points_list, error_list, 'MarkerFaceColor', method_colour, 'MarkerEdgeColor', 'none');

                    xlabel('number points [-]');
                    ylabel(error_label);
                    ylim([error_min, error_max]);
                end

                title(Tiled_Chart_Layout(m), method_label, 'FontSize', 15);
            end

            set(Axes, 'FontSize', 15);
            set(Axes, 'LineWidth', 2);

            % Save the figure
            export_fig(figure_number, [figure_name, '.fig']);
            export_fig(figure_number, [figure_name, '.png']);

            movefile([figure_name, '*'], save_folder_name); 

            %--% Relationship with uncertainty %--%
            % Average relative radius error for each cylinder with a non-empty point cloud
            fuzzy_cyl_radius_error_list = Fuzzy_Radius_Errors.Cylinders.error(~empty_bool);
            discr_cyl_radius_error_list = Discr_Radius_Errors.Cylinders.error(~empty_bool);

            % Average geometric uncertainty 
            Avg_Geom_Uncertainty_fun    = @(sigma_prop_list, sigma_radial_list) mean(sqrt(sigma_prop_list .* sigma_radial_list));
            avg_geom_uncertainty_list   = cellfun(Avg_Geom_Uncertainty_fun, sigma_prop_cell, sigma_radial_cell);

            % Plot
            figure_number   = figure_number + 1;
            Fig             = figure(figure_number);
            figure_name     = 'Radius_Error_Uncertainty';
  
            set(Fig, 'name', figure_name, 'NumberTitle', 'off');

            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85]);
            set(gcf, 'color', [1, 1, 1]);

            % Discrete relative error
            subplot(1, 2, 1);
            hold on
            grid on

            scatter(avg_geom_uncertainty_list, discr_cyl_radius_error_list, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');

            xlabel('$\bar{\sigma}$ [mm]', 'interpreter', 'latex');
            ylabel('Discrete radius rel. error [%]');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
    
            hold off

            % Fuzzy relative error
            subplot(1, 2, 2);
            hold on
            grid on

            scatter(avg_geom_uncertainty_list, fuzzy_cyl_radius_error_list, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

            xlabel('$\bar{\sigma}$ [mm]', 'interpreter', 'latex');
            ylabel('Fuzzy radius rel. error [%]');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
    
            hold off

            % Save the figure
            export_fig(figure_number, [figure_name, '.fig']);
            export_fig(figure_number, [figure_name, '.png']);

            movefile([figure_name, '*'], save_folder_name); 
        end

    %% Stem taper curve %% 
        if Compute_Radius_Errors == true
            % Taper curve data
            Stem_Taper_Diagnostics = false;
            [discr_stem_height_list, object_stem_radius_list, discrete_stem_radius_list, discrete_stem_radius_error_list, discrete_stem_radius_rel_error_list]  = Stem_Taper_Curve(Discrete_QSM, Smoothed_Model_Circle_Geometry, Parallel_Pool, Stem_Taper_Diagnostics);
            [fuzzy_stem_height_list, ~, fuzzy_stem_radius_list, fuzzy_stem_radius_error_list, fuzzy_stem_radius_rel_error_list]                                 = Stem_Taper_Curve(Fuzzy_QSM, Smoothed_Model_Circle_Geometry, Parallel_Pool, Stem_Taper_Diagnostics);
    
            % Table
            discrete_stem_mean_error        = mean(discrete_stem_radius_error_list);
            discrete_stem_mean_abs_error    = mean(abs(discrete_stem_radius_error_list));
            discrete_stem_mean_rel_error    = mean(discrete_stem_radius_rel_error_list);
    
            fuzzy_stem_mean_error           = mean(fuzzy_stem_radius_error_list);
            fuzzy_stem_mean_abs_error       = mean(abs(fuzzy_stem_radius_error_list));
            fuzzy_stem_mean_rel_error       = mean(fuzzy_stem_radius_rel_error_list);
    
            stem_error_matrix   = [discrete_stem_mean_error, discrete_stem_mean_abs_error, discrete_stem_mean_rel_error; fuzzy_stem_mean_error, fuzzy_stem_mean_abs_error, fuzzy_stem_mean_rel_error];
            table_column_labels = {'stem r mean error [m]', 'stem r mean abs error [m]', 'stem r mean rel error [%]'};
            table_file_name     = 'Stem_Taper_Curve_Errors.xls';
    
            Stem_Taper_Table = Table_Formatter(stem_error_matrix, num_decimal_digits, output_format, method_labels, table_column_labels, table_file_name, Print);

            Results_Tables.Stem_Taper = Stem_Taper_Table;
            movefile(table_file_name, save_folder_name);
        end

        % Plot
        if Plot == true && Compute_Radius_Errors == true
            % Figure
            figure_number   = figure_number + 1;
            Fig             = figure(figure_number);
            figure_name     = 'Stem_Taper_Curve';
  
            set(Fig, 'name', figure_name, 'NumberTitle', 'off');

            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85]);
            set(gcf, 'color', [1, 1, 1]);

            % First subplot showing both the taper curves
            subplot(1, 3, 1)
            hold on
            grid on

            plot(object_stem_radius_list, discr_stem_height_list, 'LineWidth', 2, 'color', 'k', 'DisplayName', 'Object');
            plot(discrete_stem_radius_list, discr_stem_height_list, 'LineWidth', 2, 'color', 'b', 'DisplayName', 'Discrete');
            plot(fuzzy_stem_radius_list, fuzzy_stem_height_list, 'LineWidth', 2, 'color', 'r', 'DisplayName', 'Fuzzy');

            xlabel('radius [m]');
            ylabel('height [m]');

            legend('show', 'location', 'northeast');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
    
            hold off

            % Second subplot for the error
            subplot(1, 3, 2)
            hold on
            grid on

            plot(discrete_stem_radius_error_list, discr_stem_height_list, 'LineWidth', 2, 'color', 'b', 'DisplayName', 'Discrete');
            plot(fuzzy_stem_radius_error_list, fuzzy_stem_height_list, 'LineWidth', 2, 'color', 'r', 'DisplayName', 'Fuzzy');

            xlabel('radius error [m]');
            ylabel('height [m]');

            legend('show', 'location', 'northeast');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
    
            hold off

            % Third subplot for the relative error
            subplot(1, 3, 3)
            hold on
            grid on

            plot(discrete_stem_radius_rel_error_list, discr_stem_height_list, 'LineWidth', 2, 'color', 'b');
            plot(fuzzy_stem_radius_rel_error_list, fuzzy_stem_height_list, 'LineWidth', 2, 'color', 'r');

            xlabel('radius rel. error [%]');
            ylabel('height [m]');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
    
            hold off

            % Save the figure
            export_fig(figure_number, [figure_name, '.fig']);
            export_fig(figure_number, [figure_name, '.png']);

            movefile([figure_name, '*'], save_folder_name); 
        end
        
    %% Saving the results %%
        % Both the fuzzy QSM and the result tables are saved in one .mat file
        save('FuzzyQSM_Outputs', 'Fuzzy_QSM', 'Results_Tables', '-v7.3');
        movefile('FuzzyQSM_Outputs.mat', save_folder_name);

end