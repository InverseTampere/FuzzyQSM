% The difference in radius between each object circle and the nearest QSM cylinder is determined
% The error is defined as QSM - object
% The Radius_Errors structure consists of Total and Branch_Order, which both contain structures for the Error and Relative_Error of the radius
% Those structures contain the data, mean and mean absolute (mean_abs) values
% Additionally the average (relative) error for each cylinder is given in the Cylinders structure

% Note that each index in Branch_Order is the branch order minus 1 and that the original object QSM's branch orders are used

function [Radius_Errors, Object_Circle_Geometry] = QSM_Radius_Error(QSM, Object_Circle_Geometry, Parallel_Pool, Plot)

    %% Inputs %%
        % Cylinder geometry
        cylinder_radius_list        = QSM.cylinder.radius;

        % Circle geometry
        circle_branch_order_list    = Object_Circle_Geometry.BranchOrder;

    %% Circle to QSM cylinder association %%
        % Association is considered in terms of nearest Euclidean distance
        QSM_Circle_Geometry     = Smoothed_Model_Circle_Cylinder_Association(QSM, Object_Circle_Geometry, Parallel_Pool);
        circle_centre_matrix    = QSM_Circle_Geometry.centre;
        circle_cyl_index_list   = QSM_Circle_Geometry.cylinder_index;
        circle_radius_list      = QSM_Circle_Geometry.radius;

    %% Radius errors %%
        % Differences between the radii
        radius_error_list       = cylinder_radius_list(circle_cyl_index_list) - circle_radius_list;
        radius_rel_error_list   = radius_error_list ./ circle_radius_list * 100;

        Radius_Errors = struct('Total', struct('error', struct('data', radius_error_list, 'mean', mean(radius_error_list), 'mean_abs', mean(abs(radius_error_list))), ...
                                               'relative_error', struct('data', radius_rel_error_list, 'mean', mean(radius_rel_error_list), 'mean_abs', mean(abs(radius_rel_error_list)))));

        % Error per cylinder
        number_cylinders        = length(cylinder_radius_list);
        cylinder_error_list     = NaN(1, number_cylinders);
        cylinder_rel_error_list = NaN(1, number_cylinders);

        for c = 1 : number_cylinders
            % Circles corresponding to this cylinder
            circle_bool = circle_cyl_index_list == c;

            % Mean error
            cylinder_error_list(c)      = mean(radius_error_list(circle_bool));
            cylinder_rel_error_list(c)  = mean(radius_rel_error_list(circle_bool));
        end

        Radius_Errors.Cylinders = struct('error', cylinder_error_list, 'relative_error', cylinder_rel_error_list);

        % Branch orders
        max_branch_order = single(max(circle_branch_order_list));       % Note conversion to single

        for branch_order = 0 : max_branch_order
            % The circles and their errors in this branch order
            branch_order_circle_bool    = circle_branch_order_list == branch_order;
            BO_radius_error_list        = radius_error_list(branch_order_circle_bool);
            BO_radius_rel_error_list    = radius_rel_error_list(branch_order_circle_bool);

            Radius_Errors.Branch_Order(branch_order + 1) = struct('error', struct('data', BO_radius_error_list, 'mean', mean(BO_radius_error_list), 'mean_abs', mean(abs(BO_radius_error_list))), ...
                                                                  'relative_error', struct('data', BO_radius_rel_error_list, 'mean', mean(BO_radius_rel_error_list), 'mean_abs', mean(abs(BO_radius_rel_error_list))));
        end
    
    %% Plot %%
        if Plot == true
            %--% Error for each circle in 3D %--%
            % Colour maps used for the errors
            number_colours  = 1e3;
            error_cmap      = cbrewer('div', 'RdBu', number_colours);
            error_cmap      = max(error_cmap, 0);
            error_cmap      = min(error_cmap, 1);

            % Plotted data sets
            error_fields    = {'error', 'relative_error'};
            error_labels    = {'radius error [m]', 'radius rel. error [%]'};
            number_errors   = length(error_fields);

            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85]);
            set(gcf, 'color', [1, 1, 1]);

            for e = 1 : number_errors
                % This error's data
                error_field = error_fields{e};
                error_label = error_labels{e};
                error_list  = Radius_Errors.Total.(error_field).data;

                max_abs_error       = max(abs(error_list));
                error_list_n        = (error_list + max_abs_error) / (2*max_abs_error);
                error_colour_ind    = round(error_list_n * (number_colours - 1)) + 1;
                error_colours       = error_cmap(error_colour_ind, :);

                % The subplot
                subplot(1, number_errors, e);
                hold on
                grid on

                sc_size = 1e1;
                scatter3(circle_centre_matrix(:, 1), circle_centre_matrix(:, 2), circle_centre_matrix(:, 3), sc_size, error_colours);

                cb = colorbar;
                shading interp
                clim(max_abs_error * [-1, 1])
                ylabel(cb, error_label);
                cb.FontSize = 15;      

                colormap(error_cmap);

                xlabel('x [m]');
                ylabel('y [m]');
                zlabel('z [m]');

                axis equal
                view(45, 45);
                
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
            end

            %--% Box plots %--%
            % Box plot settings
            box_spacing     = 1 - 0.05;
            box_width       = box_spacing - 0.05;

            figure(2)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85]);
            set(gcf, 'color', [1, 1, 1]);

            for e = 1 : number_errors
                % This error type's subplot
                error_field = error_fields{e};
                error_label = error_labels{e};

                subplot(1, number_errors, e)
                hold on
                grid on

                % Box plots
                for branch_order = 0 : max_branch_order
                    % Error data for this branch order
                    error_list = Radius_Errors.Branch_Order(branch_order + 1).(error_field).data;
        
                    % Box plot
                    boxchart(branch_order * ones(length(error_list), 1), error_list, 'BoxFaceColor', 'b', 'MarkerColor', 'b', 'BoxWidth', box_width, 'Notch', 'off');
                end
    
                % Axes
                xlim([-0.5, max_branch_order + 0.5]);
                xlabel('branch order [-]');
                ylabel(error_label);
    
                % Formatting
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
    
                hold off
            end

            % Pause message
            disp('The radius errors have been determined. The figures will close and script end when a key is pressed')
            pause()
    
            close([1, 2])
        end
end