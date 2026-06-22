% This script determines which QSM cylinder each circle corresponds to and adds it to the structure

function Smoothed_Model_Circle_Geometry = Smoothed_Model_Circle_Cylinder_Association(QSM, Smoothed_Model_Circle_Geometry, Parallel_Pool)

    %% Structure inputs %%
        % QSM cylinder geometry
        cylinder_axis_matrix        = QSM.cylinder.axis;
        cylinder_start_matrix       = QSM.cylinder.start;
        cylinder_radius_list        = QSM.cylinder.radius;
        cylinder_length_list        = QSM.cylinder.length;
        cylinder_branch_order_list  = QSM.cylinder.BranchOrder;
        cylinder_branch_index_list  = QSM.cylinder.branch;

        % Circle geometry
        circle_centre_matrix        = Smoothed_Model_Circle_Geometry.centre;
        circle_radius_list          = Smoothed_Model_Circle_Geometry.radius;
        circle_normal_vector_matrix = Smoothed_Model_Circle_Geometry.normal_vector;
        number_circles              = Smoothed_Model_Circle_Geometry.number_circles;

        % Parallel pool
        Parallel_Loop               = Parallel_Pool.Parallel_Loop;
        max_cores                   = Parallel_Pool.max_cores;
        idle_timeout                = Parallel_Pool.idle_timeout;

    %% Manual inputs %% 
        Diagnostics                 = false;            % [true, false] Diagnostics plot
   
    %% Initiate the parallel pool %%
        if Parallel_Loop == true
            % The parallel pool is started
            if isempty(max_cores)
                number_cores = feature('numcores');
            else
                number_cores = max_cores;
            end
            
            Parallel_Pool_Starter(idle_timeout, number_cores)
        else
            % Otherwise the number of cores used in parallel are set to 0 to make the parfor run as a for
            number_cores = 0;
        end

    %% Circle to cylinder association %%
        % Conversion to cell array
        [number_cylinders, num_dim] = size(cylinder_axis_matrix);
        cylinder_start_cell         = mat2cell(cylinder_start_matrix, ones(1, number_cylinders), num_dim);
        cylinder_axis_cell          = mat2cell(cylinder_axis_matrix, ones(1, number_cylinders), num_dim);
        cylinder_radius_cell        = num2cell(cylinder_radius_list);
        cylinder_length_cell        = num2cell(cylinder_length_list);

        % Finding nearest cylinders
        Cylinder_Finder_fun = @(circle_centre, circle_radius, circle_normal_vector) Cylinder_Finder(circle_centre, circle_radius, circle_normal_vector, cylinder_axis_cell, cylinder_start_cell, cylinder_radius_cell, cylinder_length_cell, Diagnostics);

        cylinder_index_list = zeros(number_circles, 1);

        DQ      = parallel.pool.DataQueue;
        tick    = 0;
        N       = number_circles;
        afterEach(DQ, @ProgressUpdate);

        parfor (c = 1 : number_circles, number_cores)
            % This circle's cylinder index
            circle_centre           = circle_centre_matrix(c, :);
            circle_radius           = circle_radius_list(c);
            circle_normal_vector    = circle_normal_vector_matrix(c, :);
            cylinder_index          = Cylinder_Finder_fun(circle_centre, circle_radius, circle_normal_vector);
            cylinder_index_list(c)  = cylinder_index;

            % Progress update
            send(DQ, c);
        end

        % Indices and branch orders added to the structure
        branch_order_list                               = cylinder_branch_order_list(cylinder_index_list);
        Smoothed_Model_Circle_Geometry.BranchOrder      = branch_order_list;

        branch_index_list                               = cylinder_branch_index_list(cylinder_index_list);
        Smoothed_Model_Circle_Geometry.branch_index     = branch_index_list;           

        Smoothed_Model_Circle_Geometry.cylinder_index   = cylinder_index_list;
                            
    %% Local functions %%
        % Function that finds the cylinder for a specific circle
        function cylinder_index = Cylinder_Finder(circle_centre, circle_radius, circle_normal_vector, cylinder_axis_cell, cylinder_start_cell, cylinder_radius_cell, cylinder_length_cell, Diagnostics)
            % Cylinder centres
            Centre_fun              = @(cylinder_start, cylinder_length, cylinder_axis) cylinder_start + cylinder_length/2*cylinder_axis;
            cylinder_centre_cell    = cellfun(Centre_fun, cylinder_start_cell, cylinder_length_cell, cylinder_axis_cell, 'UniformOutput', false);

            % Project the circle onto near cylinders
            Projection_fun              = @(cylinder_axis, cylinder_start) Point_to_Vector_Projection(circle_centre, cylinder_axis, cylinder_start);
            [~, delta_cell, omega_cell] = cellfun(Projection_fun, cylinder_axis_cell, cylinder_centre_cell, 'UniformOutput', false);

            % Instead of taking delta directly, points within the cylinder have a zero delta and otherwise to the end of the cylinder
            Cylinder_Distance_fun   = @(delta, cylinder_length) max(0, abs(delta) - cylinder_length/2);
            delta_list              = cellfun(Cylinder_Distance_fun, delta_cell, cylinder_length_cell);

            % The minimal Euclidean distance is used
            omega_list                          = vertcat(omega_cell{:});
            total_distance_list                 = sqrt(omega_list.^2 + delta_list.^2);
            [total_distance, cylinder_index]    = min(total_distance_list);

            % Diagnostics plot
            if Diagnostics == true    
                fprintf('The minimum distance is %.3g m \n', total_distance);

                figure(1)
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])    
    
                hold on
                grid on
    
                % Cylinder(s)
                number_coord    = 1e2;
                cylinder_radius = cylinder_radius_cell{cylinder_index};
                cylinder_length = cylinder_length_cell{cylinder_index};
                cylinder_centre = cylinder_centre_cell{cylinder_index};
                cylinder_axis   = cylinder_axis_cell{cylinder_index};

                [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, cylinder_length, cylinder_centre, cylinder_axis, number_coord);
                surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.50, 'LineWidth', 2, 'DisplayName', 'Cylinder');
                
                % Circle
                [circle_points, ~, ~]   = Equal_Area_Circular_Sampler(circle_radius, circle_centre, circle_normal_vector, number_coord);
                scatter3(circle_points(:, 1), circle_points(:, 2), circle_points(:, 3), 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', 'DisplayName', 'Circle');

                % Axes
                xlabel('x [m]');
                ylabel('y [m]');
                zlabel('z [m]');

                axis_scale = 2*max(circle_radius, cylinder_radius);
                xlim(circle_centre(1) + axis_scale*[-1, 1]);
                ylim(circle_centre(2) + axis_scale*[-1, 1]);
                zlim(circle_centre(3) + axis_scale*[-1, 1]);

                view(45, 45);
    
                % Formatting
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                legend('show', 'location', 'eastoutside');
    
                hold off
    
                % Pause message
                disp('The cylinder has been found for this circle. The figure closes and script continues upon a key-press.');
                pause();
    
                close(1);
            end
        end

        % Progress update
        function ProgressUpdate(~)
            % Ensures that at most every P percent is printed
            P = 1;

            % The tick is updated
            tick = tick + 1;    

            % The last tick's and current tick's progress relative to P
            progress_last   = floor((tick - 1) / N * 100 / P);
            progress        = floor(tick / N * 100 / P);

            if progress - progress_last >= 1
                fprintf('   Circle cylinder association progress: %i %% \n', progress);
            end            
        end

end