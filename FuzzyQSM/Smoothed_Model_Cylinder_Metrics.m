% This script determines the average radius, length and volume of cylinders using the circles associated with them

function Cylinder_Metrics = Smoothed_Model_Cylinder_Metrics(QSM, Object_Circle_Geometry)

    %% Structure inputs %%
        % QSM cylinder geometry
        QSM_cylinder_axis_matrix    = QSM.cylinder.axis;
        QSM_cylinder_start_matrix   = QSM.cylinder.start;
        QSM_cylinder_radius_list    = QSM.cylinder.radius;
        QSM_cylinder_length_list    = QSM.cylinder.length;
        QSM_cyl_branch_index_list   = QSM.cylinder.branch;
        QSM_cyl_branch_pos_list     = QSM.cylinder.PositionInBranch;
        QSM_cyl_branch_order_list   = QSM.cylinder.BranchOrder;

        % Circle geometry
        circle_centre_matrix        = Object_Circle_Geometry.centre;
        circle_radius_list          = Object_Circle_Geometry.radius;
        circle_normal_vector_matrix = Object_Circle_Geometry.normal_vector;
        cylinder_index_list         = Object_Circle_Geometry.cylinder_index;

    %% Manual inputs %%
        unique_margin               = 1e-3;         % [-] Margin to distinguish duplicate circles
        Diagnostics                 = false;        % [true, false]

    %% New cylinder metrics %%
        [number_cylinders, num_dim] = size(QSM_cylinder_axis_matrix);
        
        cylinder_radius_list    = zeros(number_cylinders, 1);
        cylinder_length_list    = zeros(number_cylinders, 1);
        cylinder_centre_matrix  = zeros(number_cylinders, num_dim);
        cylinder_axis_matrix    = zeros(number_cylinders, num_dim);

        for c = 1 : number_cylinders
            % Cylinder geometry
            QSM_cylinder_axis   = QSM_cylinder_axis_matrix(c, :);
            QSM_cylinder_start  = QSM_cylinder_start_matrix(c, :);
            QSM_cylinder_radius = QSM_cylinder_radius_list(c);
            QSM_cylinder_length = QSM_cylinder_length_list(c);
            QSM_cylinder_centre = QSM_cylinder_start + QSM_cylinder_length/2 * QSM_cylinder_axis;

            % Circles associated with this cylinder
            circle_bool = cylinder_index_list == c;

            cyl_circle_centre_matrix        = circle_centre_matrix(circle_bool, :);
            cyl_circle_radius_list          = circle_radius_list(circle_bool);
            cyl_circle_normal_vector_matrix = circle_normal_vector_matrix(circle_bool, :);

            % Check that they are unique
            [cyl_circle_centre_matrix, number_circles, duplicate_bool]  = Unique_Margin(cyl_circle_centre_matrix, unique_margin);
            [cyl_circle_radius_list, cyl_circle_normal_vector_matrix]   = deal(cyl_circle_radius_list(~duplicate_bool), cyl_circle_normal_vector_matrix(~duplicate_bool, :));

            % If there are no associated circles, the original properties are kept
            if number_circles == 0
                cylinder_axis   = QSM_cylinder_axis;
                cylinder_centre = QSM_cylinder_centre;
                cylinder_radius = QSM_cylinder_radius;
                cylinder_length = QSM_cylinder_length;

            % If there is only one circle, only the radius is altered
            elseif number_circles == 1
                cylinder_axis   = QSM_cylinder_axis;
                cylinder_centre = QSM_cylinder_centre;
                cylinder_length = QSM_cylinder_length;
                cylinder_radius = cyl_circle_radius_list;

            % Otherwise, the circles are used to determine the average geometry
            else
                % Projection onto the QSM cylinder axis
                [~, delta_list, ~] = Point_to_Vector_Projection(cyl_circle_centre_matrix, QSM_cylinder_axis, QSM_cylinder_start);

                % The contribution of each circle is weighted by the distance to the circles next to it
                [delta_list, order]                 = sort(delta_list);
                distance_list                       = diff(delta_list);

                weight_list                         = zeros(number_circles, 1);
                weight_list(1)                      = distance_list(1);
                weight_list(number_circles)         = distance_list(number_circles - 1);
                weight_list(2 : number_circles - 1) = (distance_list(1 : number_circles - 2) + distance_list(2 : number_circles - 1)) / 2;
                weight_list                         = weight_list / sum(weight_list);

                cylinder_centre = sum(weight_list .* cyl_circle_centre_matrix(order, :), 1);
                cylinder_radius = sum(weight_list .* cyl_circle_radius_list(order));
                cylinder_axis   = sum(weight_list .* cyl_circle_normal_vector_matrix(order, :), 1);
                cylinder_axis   = cylinder_axis / norm(cylinder_axis);

                % The length follows from the two most extreme circles
                distance_matrix = pdist2(cyl_circle_centre_matrix, cyl_circle_centre_matrix);
                cylinder_length = max(distance_matrix, [], 'all');
            end

            % Adding to the matrix
            cylinder_axis_matrix(c, :)      = cylinder_axis;
            cylinder_centre_matrix(c, :)    = cylinder_centre;
            cylinder_radius_list(c)         = cylinder_radius;
            cylinder_length_list(c)         = cylinder_length;

            % Diagnostics plot
            if Diagnostics == true
                % Number of coordinates for the cylinders and circles
                number_coord = 1e2;

                figure(1)
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])    
    
                hold on
                grid on

                % Original cylinder
                [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(QSM_cylinder_radius, QSM_cylinder_length, QSM_cylinder_centre, QSM_cylinder_axis, number_coord);
                surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.50, 'LineWidth', 2, 'DisplayName', 'Original cylinder');

                % New cylinder
                [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, cylinder_length, cylinder_centre, cylinder_axis, number_coord);
                surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.50, 'LineWidth', 2, 'DisplayName', 'New cylinder');

                % Circles
                for i = 1 : number_circles
                    circle_radius           = cyl_circle_radius_list(i);
                    circle_centre           = cyl_circle_centre_matrix(i, :);
                    circle_normal_vector    = cyl_circle_normal_vector_matrix(i ,:);
                    [circle_points, ~, ~]   = Equal_Area_Circular_Sampler(circle_radius, circle_centre, circle_normal_vector, number_coord);
                    sc_circ = scatter3(circle_points(:, 1), circle_points(:, 2), circle_points(:, 3), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'DisplayName', 'Circle');

                    if i > 1
                        sc_circ.HandleVisibility = 'Off';
                    end
                end

                % Axes
                xlabel('x [m]');
                ylabel('y [m]');
                zlabel('z [m]');
                axis equal

                view(45, 45);
    
                % Formatting
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                legend('show', 'location', 'eastoutside');
    
                hold off
    
                % Pause message
                disp('The new cylinder metrics have been determined. The figure closes and script continues upon a key-press.');
                pause();
    
                close(1);
            end
        end

        % Structure for the new cylinder metrics
        cylinder_start_matrix   = cylinder_centre_matrix - 1/2*cylinder_length_list .* cylinder_axis_matrix;
        Cylinder_Metrics        = struct('radius', cylinder_radius_list, 'start', cylinder_start_matrix, 'centre', cylinder_centre_matrix, 'axis', cylinder_axis_matrix, 'length', cylinder_length_list, 'branch', QSM_cyl_branch_index_list, 'PositionInBranch', QSM_cyl_branch_pos_list, 'BranchOrder', QSM_cyl_branch_order_list);
        
        % Ensuring the cylinder interfaces are reasonable
        Cylinder_Metrics        = Branch_Interface_Fitting(Cylinder_Metrics);

end