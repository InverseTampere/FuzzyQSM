% This script takes an initial cylinder and its point cloud and fits an infinite fuzzy cylinder to it

function [Fuzzy_Inf_Cylinder, Point_Cloud_Distributions] = Fuzzy_QSM_Infinite_Cylinder_Fitting(Init_Cylinder, Point_Cloud_Coord, Fitting_Parameters, Scanner_Parameters, Scanning_Parameters)

    %% Structure inputs %%
        % Initial cylinder
        init_cylinder_start     = double(Init_Cylinder.start);
        init_cylinder_axis      = double(Init_Cylinder.axis);
        init_cylinder_radius    = double(Init_Cylinder.radius);
        init_cylinder_length    = double(Init_Cylinder.length);

        % Point cloud
        point_cloud_cell        = Point_Cloud_Coord.point_cloud_cell;

        % Fitting Parameters
        bounds_margin           = Fitting_Parameters.bounds_margin;

        % Scanner parameters
        range_bias              = Scanner_Parameters.range_bias;

        % Scanning parameters
        Scanner_loc_cell        = Scanning_Parameters.Scanner_loc_cell;
        number_scanners         = Scanning_Parameters.number_scanners;

    %% Manual inputs %%
        % If desired, intermediate results can be checked
        Transform_Diagnostics   = false;        % [true, false]
        Distr_Diagnostics       = false;        % [true, false]

        % Optimiser inputs
        max_iterations          = 1e2;          % [-] Maximum number of iterations by the optimiser. If it is reached, a warning is displayed
        function_tolerance      = 1e-4;         % [-] Minimum size for the step in the objective function. Note that it is normalised by the initial value
        step_tolerance          = 1e-4;         % [-] Minimum size for the step in the design parameters. Note that they are normalised between [0, 1]
        FD_step_size            = 1e-4;         % [-] Step size taken to evaluate the finite difference. Note that sqrt(eps) is the default
       
        % Displayed results
        Print                   = false;         % [true, false]
        Plot                    = false;         % [true, false]

    %% Fuzzy cloud %%
        % Uncertainty based on the initial cylinder geometry
        [sigma_radial_cell, sigma_prop_cell, ~, ~, ~] = Cylindrical_Object_Uncertainty(init_cylinder_start, init_cylinder_axis, Scanner_loc_cell, Scanner_Parameters, Point_Cloud_Coord);

        % The 3D distributions
        alpha                       = 0.99;          % The alpha value is irrelevant
        Point_Cloud_Distributions   = Point_Cloud_Multivariate_Normal_Generation(alpha, Point_Cloud_Coord, sigma_radial_cell, sigma_prop_cell, Scanner_loc_cell, range_bias, Distr_Diagnostics);

    %% Translation and rotation for robustness %%
        % To increase robustness, everything is centered on the centroid of the given point cloud
        [Point_Cloud_Coord_c, Point_Cloud_Distributions_c, ~, point_cloud_centroid, Scanning_Parameters_c, ~, ~] = Point_Cloud_Centering(Point_Cloud_Coord, Point_Cloud_Distributions, Scanning_Parameters);
        
        % Additionally, the data is rotated s.t. the elevation angle of the cylinder axis is pi/4 to avoid gimbal lock and such that the initial elevation angle is in the middle of the optimisation range
        num_dim                         = length(init_cylinder_axis);
        origin                          = zeros(1, num_dim);
        [~, cylinder_vector_basis, ~]   = Vector_Based_Rotation(origin, init_cylinder_axis, origin);
        
        pitch                           = pi/4;
        pitched_vector_basis            = [cos(pitch), 0, sin(pitch); 0, 1, 0; -sin(pitch), 0, cos(pitch)];

        gimbal_rotation_matrix          = pitched_vector_basis * cylinder_vector_basis;

        [Point_Cloud_Coord_r, Point_Cloud_Distributions_r, ~, Scanning_Parameters_r, Scanner_loc_cell_r, ~] = Point_Cloud_Rotation(gimbal_rotation_matrix, Point_Cloud_Coord_c, Point_Cloud_Distributions_c, Scanning_Parameters_c);

        % The initial geometry must be translated and rotated as well
        init_cylinder_axis_r                    = (gimbal_rotation_matrix * init_cylinder_axis')';
        init_cylinder_start_c                   = init_cylinder_start - point_cloud_centroid;
        init_cylinder_start_r                   = (gimbal_rotation_matrix * init_cylinder_start_c')';
        init_cylinder_centre_r                  = init_cylinder_start_r + init_cylinder_length/2 * init_cylinder_axis_r;
        [init_circle_centre_matrix_r, ~, ~, ~]  = Circle_Cylinder_Centre_Conversion([], [], init_cylinder_centre_r, init_cylinder_axis_r, origin, Scanner_loc_cell_r);
        init_circle_centre_r                    = init_circle_centre_matrix_r(1, :);

        % Diagnostic plot to check the transformation
        if Transform_Diagnostics == true
            % Scanner colour map                    
            scanner_colours = cbrewer('qual', 'Set2', max(number_scanners, 3));     % The colours used for the scanners
            scanner_colours = max(scanner_colours, 0);
            scanner_colours = min(scanner_colours, 1);
            
            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            %--% Original coordinate frame %--%
            subplot(1, 2, 1)
            hold on
            grid on
                
            % The initial cylinder surface
            number_coord                                                                = 1e3;     
            init_cylinder_centre                                                        = init_cylinder_start + init_cylinder_length/2 * init_cylinder_axis;
            [init_cylinder_coord_x, init_cylinder_coord_y, init_cylinder_coord_z, ~]    = Cylinder_Surface_Generator(init_cylinder_radius, init_cylinder_length, init_cylinder_centre, init_cylinder_axis, number_coord);
            surf(init_cylinder_coord_x, init_cylinder_coord_y, init_cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'm', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Initial cylinder');
  
            % The point cloud        
            for s = 1 : number_scanners
                % This scanner's point cloud
                scanner_colour          = scanner_colours(s, :);
                point_cloud_matrix    = point_cloud_cell{s};

                scatter3(point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 2), 'MarkerFaceColor', scanner_colour, 'MarkerEdgeColor', 'none', 'DisplayName', 'Point cloud');
            end

            % Axes
            xlabel('x [m]')
            ylabel('y [m]')
            zlabel('z [m]')
            axis equal
            
            % Viewing angle
            view(45, 45)

            % Legend
            legend('show', 'location', 'northoutside');

            % Font size
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            %--% Transformed coordinate frame %--%
            subplot(1, 2, 2)
            hold on
            grid on
                
            % The initial cylinder surface
            [init_cylinder_coord_a, init_cylinder_coord_b, init_cylinder_coord_c, ~] = Cylinder_Surface_Generator(init_cylinder_radius, init_cylinder_length, init_cylinder_centre_r, init_cylinder_axis_r, number_coord);
            surf(init_cylinder_coord_a, init_cylinder_coord_b, init_cylinder_coord_c, 'EdgeColor', 'none', 'FaceColor', 'm', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Initial cylinder');
  
            % The point cloud
            cyl_point_cloud_cell = Point_Cloud_Coord_r.point_cloud_cell;

            for s = 1 : number_scanners
                % This scanner's point cloud
                scanner_colour          = scanner_colours(s, :);
                point_cloud_matrix    = cyl_point_cloud_cell{s};
                
                scatter3(point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 2), 'MarkerFaceColor', scanner_colour, 'MarkerEdgeColor', 'none', 'DisplayName', 'Point cloud');
            end

            % Axes
            xlabel('a [m]')
            ylabel('b [m]')
            zlabel('c [m]')
            axis equal
            
            % Viewing angle
            view(45, 45)

            % Legend
            legend('show', 'location', 'northoutside');

            % Font size
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off          
        end

   %% Solution bounds %%
        % The geometry parameters' indices (as the optimiser's inputs must be vector format)
        centre_ind                  = 1 : 2;
        azim_ind                    = max(centre_ind) + 1;
        elev_ind                    = azim_ind + 1;
        radius_ind                  = elev_ind + 1;
        num_optim_geo_parameters    = radius_ind;

        Geometry_Indices            = struct('centre', centre_ind, 'azim_angle', azim_ind, 'elev_angle', elev_ind, 'radius', radius_ind);

        % Finally, the set of initial centered geometry parameters
        inf_cylinder_geometry_init                                                  = zeros(1, num_optim_geo_parameters);
        [azim_angle_init_r, elev_angle_init_r]                                      = deal(0, pitch);
        inf_cylinder_geometry_init([centre_ind, azim_ind, elev_ind, radius_ind])    = [init_circle_centre_r, azim_angle_init_r, elev_angle_init_r, init_cylinder_radius];

        % Reasonable bounds for the optimisation process are based on the least-squares geometry s.t. the normalised variables can be directly compared
        circle_centre_LB = init_circle_centre_r - bounds_margin * init_cylinder_radius;
        circle_centre_UB = init_circle_centre_r + bounds_margin * init_cylinder_radius;
        
        radius_LB = init_cylinder_radius / (1 + bounds_margin);
        radius_UB = init_cylinder_radius * (1 + bounds_margin);

        % The direction vectors can deviate by the radius times the margin, at the cylinder length's distance
        % Note that this deviation only results in allowed changes to the elevation angle. The azimuth angle remains unbounded
        elev_deviation      = atan(init_cylinder_radius * bounds_margin / init_cylinder_length);

        elev_LB             = max(-pi/2, elev_angle_init_r - elev_deviation);
        elev_UB             = min(pi, elev_angle_init_r + elev_deviation);
        [azim_LB, azim_UB]  = deal(-3*pi, 3*pi);

        % The combined sets of bounds
        geometry_LB                                                 = zeros(1, num_optim_geo_parameters);
        geometry_LB([centre_ind, azim_ind, elev_ind, radius_ind])   = [circle_centre_LB, azim_LB, elev_LB, radius_LB];

        geometry_UB                                                 = zeros(1, num_optim_geo_parameters);
        geometry_UB([centre_ind, azim_ind, elev_ind, radius_ind])   = [circle_centre_UB, azim_UB, elev_UB, radius_UB];

        % As a result, the optimisation bounds are between 0 and 1
        Optim_LB = zeros(1, num_optim_geo_parameters);
        Optim_UB = ones(1, num_optim_geo_parameters);
        
    %% Infinite cylinder %%
        %--% Optimisation %--%
        % Centered geometry and a centered point cloud are used with the initial objective value for normalisation
        inf_cylinder_geometry_init_n    = (inf_cylinder_geometry_init - geometry_LB) ./ (geometry_UB - geometry_LB);
        inf_cyl_objective_value_init    = Infinite_Cylinder_Objective_Function(inf_cylinder_geometry_init_n, Geometry_Indices, geometry_LB, geometry_UB, 1, Point_Cloud_Distributions_r, Scanning_Parameters_r);

        Inf_Cyl_Objective_Function_fun  = @(inf_cylinder_geometry_n) Infinite_Cylinder_Objective_Function(inf_cylinder_geometry_n, Geometry_Indices, geometry_LB, geometry_UB, inf_cyl_objective_value_init, Point_Cloud_Distributions_r, Scanning_Parameters_r);

        % Certain information is saved during optimisation
        [Infinite_cylinder_geometry_steps_n, Gradient_steps_n, Objective_value_steps] = deal([]);

        % Optimisation with the interior-point algorithm
        Options = optimoptions('fmincon', 'Display', 'off', 'OutputFcn', @Optim_Path, 'Algorithm', 'interior-point', 'SubProblemAlgorithm', 'cg', 'EnableFeasibilityMode', true, 'MaxIterations', max_iterations, 'FunctionTolerance', function_tolerance, 'StepTolerance', step_tolerance, 'FiniteDifferenceStepSize', FD_step_size);

        [inf_cylinder_geometry_n, Objective_value, ~, Output, ~, ~, Optimum_Hessian_matrix_n] = fmincon(Inf_Cyl_Objective_Function_fun, inf_cylinder_geometry_init_n, [], [], [], [], Optim_LB, Optim_UB, [], Options);

        %--% Checks %--%
        % Check if the optimum is worse than the initial geometry, in which case sqp is used as the optimisation algorithm
        if Objective_value > 1
            % The optimisation data is cleared
            [Infinite_cylinder_geometry_steps_n, Gradient_steps_n, Objective_value_steps] = deal([]);

            % Updated optimisation method
            Options.Algorithm = 'sqp';
            
            % New optimisation run
            [inf_cylinder_geometry_n, Objective_value, ~, Output, ~, ~, Optimum_Hessian_matrix_n] = fmincon(Inf_Cyl_Objective_Function_fun, inf_cylinder_geometry_init_n, [], [], [], [], Optim_LB, Optim_UB, [], Options);    % Optimisation with sqp
        end
        
        % Check if the maximum number of iterations was reached
        number_optimiser_steps = Output.iterations;

        if number_optimiser_steps == max_iterations
            warning('optimiser:max_iter', 'The maximum number of iterations for the optimiser was reached. Check whether this should be increased or whether the objective function is poor.')
        end

        % Check if any of the bounds are reached, which may indicate that the optimum lies outside the optimisation space
        [centre_n, azim_n, elev_n, radius_n] = deal(inf_cylinder_geometry_n(centre_ind), inf_cylinder_geometry_n(azim_ind), inf_cylinder_geometry_n(elev_ind), inf_cylinder_geometry_n(radius_ind));

        if max(centre_n) > 1 - 2*step_tolerance || min(centre_n) < 2*step_tolerance                                                                             % Note that the finite step size means the bounds are unlikely to be hit exactly
            warning('bounds:centre', 'The bounds for optimisation have restricted the optimal centre location. Try increasing the bounds margin.');
        end
        if radius_n > 1 - 2*step_tolerance || radius_n < 2*step_tolerance
            warning('bounds:radius', 'The bounds for optimisation have restricted the optimal radius. Try increasing the bounds margin.');
        end
        if (azim_n < 2*step_tolerance && azim_LB > -pi) || (azim_n > 1 - 2*step_tolerance && azim_UB < pi) || ...
           (elev_n < 2*step_tolerance && elev_LB > 0) || (elev_n > 1 - 2*step_tolerance && elev_UB < pi/2)                                                      % The bounds are allowed to be hit if the angles are the spherical angle bounds
            warning('bounds:axis', 'The bounds for optimisation have restricted the optimal cylinder axis direction. Try increasing the bounds margin.');
        end

        % Check if the objective function is convex at the solution (i.e. if all eigenvalues of the Hessian are positive)
        % Note that as the optimising algorithm is quasi-Newton, convexity is not guaranteed
        Hessian_n_eigenvalues       = eig(Optimum_Hessian_matrix_n);
        Hessian_n_eigenvalues_sign  = sign(Hessian_n_eigenvalues);

        if sum(Hessian_n_eigenvalues_sign) < num_optim_geo_parameters
            warning('optimiser:convexity', 'The objective function is not convex at the found solution. The lowest Hessian eigenvalue is %.3g \n', min(Hessian_n_eigenvalues));
        end

        % Cross-sectional geometry
        [cylinder_radius, ~, cylinder_direction, ~, cylinder_centre] = Geometry_Retrieval(inf_cylinder_geometry_n, geometry_LB, geometry_UB, gimbal_rotation_matrix, point_cloud_centroid, Scanner_loc_cell_r, Geometry_Indices);

        % Structure containing the infinite cylinder geometry
        Fuzzy_Inf_Cylinder = struct('radius', cylinder_radius, 'axis', cylinder_direction, 'centre', cylinder_centre);

    %% Fitting information %%
        % Printed statement
        if Print == true
            % The objective value is converted into the expected Mahalanobis distance
            Expected_Mahal_distance = Objective_value * inf_cyl_objective_value_init;

            fprintf('The infinite cylinder has been fitted with O =  %.4g, E[M^2] = %.4g \n', Objective_value, Expected_Mahal_distance);
            fprintf('   Radius:     %.3g m \n',             cylinder_radius);
            fprintf('   Centre:     [%.3g, %.3g, %.3g] \n', cylinder_centre);
            fprintf('   Direction:  [%.3g, %.3g, %.3g] \n', cylinder_direction);
            fprintf('It took a total of %i steps \n', number_optimiser_steps);
        end      
        
        % Plots showing the point cloud, fitted cylinder and optimisation history
        if Plot == true
            %--% Geometry overview %--%
            % The number of coordinates in each distribution (if shown) and the cylinder
            number_coord        = 1e2;     
            max_distributions   = 1e2;
        
            % Colours for the point cloud of each scanner
            scanner_colours = cbrewer('qual', 'Set2', max(number_scanners, 3));     % The colours used for the scanners
            scanner_colours = max(scanner_colours, 0);
            scanner_colours = min(scanner_colours, 1);
            
            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            hold on
            grid on
                
            % The initial cylinder surface
            init_cylinder_centre                                                        = init_cylinder_start + init_cylinder_length/2 * init_cylinder_axis;
            [init_cylinder_coord_x, init_cylinder_coord_y, init_cylinder_coord_z, ~]    = Cylinder_Surface_Generator(init_cylinder_radius, init_cylinder_length, init_cylinder_centre, init_cylinder_axis, number_coord);
            surf(init_cylinder_coord_x, init_cylinder_coord_y, init_cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'm', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Initial cylinder');
            
            % The optimised cylinder surface
            [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, init_cylinder_length, cylinder_centre, cylinder_direction, number_coord);
            surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Fuzzy cylinder');

            % The point cloud and its uncertainty
            distribution_mu_cell    = Point_Cloud_Distributions.distribution_mu_cell;
            distribution_radii_cell = Point_Cloud_Distributions.distribution_radii_cell;
            distribution_axes_cell  = Point_Cloud_Distributions.distribution_axes_cell;
            
            number_distributions_list   = Point_Cloud_Distributions.number_distributions_list;
            num_cumulative_distr_list   = [0, cumsum(number_distributions_list)];
            
            for s = 1 : number_scanners
                %--% Scanner vectors %--%
                % Distance from scanner to cylinder
                scanner_loc     = Scanner_loc_cell{s};
                scanner_vector  = scanner_loc - cylinder_centre;
                range           = sqrt(sum(scanner_vector.^2));
                
                % The vector
                scanner_vector_scaled   = init_cylinder_length * scanner_vector / norm(scanner_vector);
                scanner_string          = sprintf('Scanner %g, R = %.3g m', s, range);
                scanner_colour          = scanner_colours(s, :);
                plot3(cylinder_centre(1) + [0, scanner_vector_scaled(1)], cylinder_centre(2) + [0, scanner_vector_scaled(2)], cylinder_centre(3) + [0, scanner_vector_scaled(3)], 'LineWidth', 2, 'color', scanner_colour, 'DisplayName', scanner_string);
                
                %--% Scanner's point cloud %--%
                number_distributions_s = number_distributions_list(s);

                start_ind   = num_cumulative_distr_list(s) + 1;
                end_ind     = num_cumulative_distr_list(s + 1);

                if number_distributions_s > max_distributions
                    % If there are too many distributions, only mu is shown
                    distribution_mu_matrix = vertcat(distribution_mu_cell{start_ind : end_ind});

                    scatter3(distribution_mu_matrix(:, 1), distribution_mu_matrix(:, 2), distribution_mu_matrix(:, 3), 'MarkerFaceColor', scanner_colour, 'MarkerEdgeColor', 'none', 'DisplayName', sprintf('%s', '\mu'));
                else                    
                    for d = start_ind : end_ind
                        % Distribution properties
                        distribution_mu     = distribution_mu_cell{d};
                        distribution_radii  = distribution_radii_cell{d};
                        distribution_axes   = distribution_axes_cell{d};
    
                        % Surface
                        [ellipsoid_coord_matrix, number_coord] = Ellipsoid_Coordinate_Generator(distribution_mu, distribution_radii, distribution_axes, number_coord);
    
                        x_ellipsoid = reshape(ellipsoid_coord_matrix(:, 1), sqrt(number_coord) * [1, 1]);
                        y_ellipsoid = reshape(ellipsoid_coord_matrix(:, 2), sqrt(number_coord) * [1, 1]);
                        z_ellipsoid = reshape(ellipsoid_coord_matrix(:, 3), sqrt(number_coord) * [1, 1]);
    
                        surf_el = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', scanner_colour, 'FaceAlpha', 0.10, 'DisplayName', sprintf('Point cloud %g uncertainty, \\alpha = %.3g', s, alpha));
    
                        % Mu
                        sc_mu = scatter3(distribution_mu(1), distribution_mu(2), distribution_mu(3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', scanner_colour, 'DisplayName', '\mu');
    
                        if d > start_ind
                            sc_mu.HandleVisibility      = 'Off';
                            surf_el.HandleVisibility    = 'Off';
                        end
                    end 
                end
            end
            
            % Aspect ratio
            axis equal

            % Axes
            xlabel('x [m]')
            ylabel('y [m]')
            zlabel('z [m]')
            
            % Viewing angle
            view(45, 45)

            % Legend
            legend('show', 'location', 'eastoutside');

            % Font size
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off
                        
            %--% Optimisation parameter history %--%
            % Step sizes (based on the normalised geometry)
            Geometry_steps_n_delta_matrix       = diff(Infinite_cylinder_geometry_steps_n, 1, 1);
            Optimiser_step_size_list            = [0; sqrt(sum(Geometry_steps_n_delta_matrix.^2, 2))];              % A 0 value is appended for the first step
            
            optim_parameter_labels  = {'Q', '\Delta'};
            optim_parameter_units   = {'-', '-'};
            num_optim_parameters    = length(optim_parameter_labels);
            
            optim_parameter_colours = cbrewer('qual', 'Set2', max(num_optim_parameters, 3));
            optim_parameter_colours = max(optim_parameter_colours, 0);
            optim_parameter_colours = min(optim_parameter_colours, 1);
            
            optim_parameter_history_matrix = [Objective_value_steps, Optimiser_step_size_list];
            
            figure(2)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])   
            
            for q = 1 : num_optim_parameters
                % This parameter's data
                optim_parameter_label   = optim_parameter_labels{q};
                optim_parameter_unit    = optim_parameter_units{q};
                optim_parameter_colour  = optim_parameter_colours(q, :);
                optim_parameter_history = optim_parameter_history_matrix(:, q);
                
                % The subplot
                subplot(1, num_optim_parameters, q)
                hold on
                grid on

                plot(1 : length(optim_parameter_history), optim_parameter_history, 'LineWidth', 2, 'Color', optim_parameter_colour);

                % Axes
                xlim([1, length(optim_parameter_history)]);
                xlabel('Optim. step [-]');
                ylabel(sprintf('%s [%s]', optim_parameter_label, optim_parameter_unit));

                % Font
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off   
            end

            %--% Pause message %--%
            disp('The script will continue and the plots will close when a key is pressed');
            pause();

            close all;             
        end

    %% Local functions %%
        % Infinite cylinder fitting objective function
        function Objective_value = Infinite_Cylinder_Objective_Function(inf_cylinder_geometry_n, Geometry_Indices, geometry_LB, geometry_UB, Objective_value_init, Point_Cloud_Distributions_c, Scanning_Parameters_c)

            %% Structure inputs %%
                % Geometry indices
                centre_indices              = Geometry_Indices.centre;
                radius_index                = Geometry_Indices.radius;
                azim_index                  = Geometry_Indices.azim_angle;
                elev_index                  = Geometry_Indices.elev_angle;

                % Distributions
                number_distributions        = Point_Cloud_Distributions_c.number_distributions;
                distribution_mu_cell_c      = Point_Cloud_Distributions_c.distribution_mu_cell;
                distribution_axes_cell      = Point_Cloud_Distributions_c.distribution_axes_cell;
                distribution_sigmae_cell    = Point_Cloud_Distributions_c.distribution_sigmae_cell;
        
                % Scanner locations
                Scanner_loc_cell_c          = Scanning_Parameters_c.Scanner_loc_cell;

            %% Manual inputs %%
                % Fitting parameters
                distance_moment             = 2;            % [-] E[M^d]

                % It is recommend that outputs are set to false        
                Mahal_Diagnostics           = false;        % [true, false] Diagnostics plot for the expected Mahalanobis distance
                Objective_Diagnostics       = false;        % [true, false] Shows the objective values for each distribution

            %% Geometry parameters (optimisation variables) %%
                % Un-normalising the geometry vector
                inf_cylinder_geometry = inf_cylinder_geometry_n .* (geometry_UB - geometry_LB) + geometry_LB;
        
                % Circular cross section
                circle_centre   = inf_cylinder_geometry(centre_indices);               % Expressed in the cylinder-aligned coordinate frame (i.e. projection on the plane orthogonal to the cylinder axis crossing the origin)
                circle_radius   = inf_cylinder_geometry(radius_index);    
                
                % Cylinder direction from the spherical angles
                azim_angle              = inf_cylinder_geometry(azim_index);
                azim_angle              = azim_angle - sign(azim_angle)*2*pi;
        
                elev_angle              = inf_cylinder_geometry(elev_index);
                
                if elev_angle < 0
                    elev_angle          = abs(elev_angle);
                elseif elev_angle > pi/2
                    elev_angle          = pi - elev_angle;
                end
        
                [cylinder_dir, ~, ~]    = Vector_Spherical_Angle_Conversion([], azim_angle, elev_angle);
        
                % The circle centre is converted to the cylinder centre
                num_dim = length(cylinder_dir);
                origin  = zeros(1, num_dim);                                                            % The origin is used as the centroid as the point cloud is already centered
                
                circle_height                   = 0;                                                    % Note that the height above the cross-section is irrelevant when the problem is projected onto it
                [~, ~, ~, cylinder_centre_c]    = Circle_Cylinder_Centre_Conversion(circle_centre, circle_height, [], cylinder_dir, origin, Scanner_loc_cell_c);
    
            %% The objective value %%
                % The expected Mahalanobis distances are the objective value for each distribution
                [~, distr_objective_list] = Expected_Mahal_Distance_Cylinder_Line_Approx(cylinder_centre_c, circle_radius, cylinder_dir, distance_moment, Point_Cloud_Distributions_c, Scanner_loc_cell_c, Mahal_Diagnostics);

                % The final weighted objective value
                Objective_value_list    = distr_objective_list / Objective_value_init;
                Objective_value         = mean(Objective_value_list);
    
            %% Diagnostics %%
                % For diagnostic purposes the objective value and geometry are printed
                if Objective_Diagnostics == true
                    % Printed messages
                    disp('---');
                    fprintf('O = %.3g \n', Objective_value);
                    fprintf('c = [%.3g, %.3g, %.3g] m \n', cylinder_centre_c);
                    fprintf('v = [%.3g, %.3g, %.3g] m \n', cylinder_dir);
                    fprintf('r = %.3g m \n', circle_radius);       
        
                    % Objective colour map
                    number_colours  = 1e3;
                    objective_cmap = cbrewer('seq', 'Blues', number_colours);
                    objective_cmap = max(objective_cmap, 0);
                    objective_cmap = min(objective_cmap, 1);
        
                    min_objective       = min(distr_objective_list_w);
                    max_objective       = max([distr_objective_list_w; min_objective + 1e-16]);                             % In case the values are all identical
                    objective_list_n    = (distr_objective_list_w - min_objective) / (max_objective - min_objective);
                    cmap_ind            = round((number_colours - 1) * objective_list_n) + 1;
        
                    % Number of coordinates per cylinder/distribution
                    number_coord = 1e2;
        
                    % Plot
                    figure(1)
                    % Set the size and white background color
                    set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                    set(gcf, 'color', [1, 1, 1])    
        
                    % First axes for the cylinder and distribution
                    ax1 = axes;
                    hold on
                    grid on
        
                    % Cylinder
                    distribution_mu_matrix_c    = vertcat(distribution_mu_cell_c{:});
                    distance_list               = sqrt(sum((distribution_mu_matrix_c - cylinder_centre_c).^2, 2));
                    cylinder_length_margin      = 2*max(distance_list);
        
                    [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(circle_radius, cylinder_length_margin, cylinder_centre_c, cylinder_dir, number_coord);
                    surf_cyl_name   = 'Cylinder';
                    surf_cyl        = surf(ax1, cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', surf_cyl_name);
        
                    % Distributions
                    for i = 1 : number_distributions
                        % Distribution properties
                        distribution_mu     = distribution_mu_cell_c{i};
                        distribution_axes   = distribution_axes_cell{i};
                        distribution_sigmae = distribution_sigmae_cell{i};
        
                        % Coordinates
                        [distr_coord_matrix, number_coord] = Ellipsoid_Coordinate_Generator(distribution_mu, distribution_sigmae, distribution_axes, number_coord);
            
                        x_distribution  = reshape(distr_coord_matrix(:, 1), sqrt(number_coord) * [1, 1]);
                        y_distribution  = reshape(distr_coord_matrix(:, 2), sqrt(number_coord) * [1, 1]);
                        z_distribution  = reshape(distr_coord_matrix(:, 3), sqrt(number_coord) * [1, 1]);
            
                        % Distribution at 1 sigma
                        distr_surf_name = sprintf('Distribution, 1%s', '\sigma');
                        distr_surf      = surf(ax1, x_distribution, y_distribution, z_distribution, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.15, 'DisplayName', distr_surf_name);
                    end
        
                    % Second axes for the points and their objective values
                    ax2 = axes;
                    hold on
                    grid on
        
                    sc_points_name = 'Point cloud';
                    sc_points      = scatter3(ax2, distribution_mu_matrix_c(:, 1), distribution_mu_matrix_c(:, 2), distribution_mu_matrix_c(:, 3), 200, objective_cmap(cmap_ind, :), 'Marker', '.', 'DisplayName', sc_points_name);
        
                    % Colour bar
                    colormap(ax2, objective_cmap);
        
                    cb = colorbar(ax2);
                    shading interp
                    clim([min_objective, max_objective])
                    ylabel(cb, 'O [-]');
                    cb.FontSize = 15;        
            
                    % Axes
                    xlabel(ax1, 'x [m]');
                    ylabel(ax1, 'y [m]');
                    zlabel(ax1, 'z [m]');
                
                    ax1.DataAspectRatio = [1, 1, 1];
                    ax2.DataAspectRatio = [1, 1, 1];
                        
                    % Legend
                    legend(ax1, [surf_cyl, distr_surf, sc_points], {surf_cyl_name, distr_surf_name, sc_points_name}, 'location', 'northoutside');
            
                    set([ax1, ax2], 'FontSize', 15);
                    set([ax1, ax2], 'LineWidth', 2);
                
                    % Change the size, link the axes and make the second set invisible
                    axis_dimensions     = [ax1.Position; ax2.Position];
                    axis_starts         = max(axis_dimensions(:, 1:2), [], 1);
                    axis_sizes          = min(axis_dimensions(:, 3:4), [], 1);    
                    ax1.Position        = [axis_starts, axis_sizes];
                    ax2.Position        = [axis_starts, axis_sizes];
                    ax2.Visible         = 'off';
            
                    linkprop([ax1, ax2], {'View', 'XLim', 'YLim', 'ZLim'});
        
                    view(45, 45);
        
                    % Pause message
                    disp('The objective values have been determined. The figure will close and script end upon a key-press.');
                    pause();
            
                    close(1);
                end
        end
    
        % Optimisation path function
        function stop = Optim_Path(inf_cylinder_geometry_n, Optimisation_Values, ~)
            stop = false;

            % Concatenate objective value, geometry and various other terms computed during each iteration
            Infinite_cylinder_geometry_steps_n  = [Infinite_cylinder_geometry_steps_n; inf_cylinder_geometry_n];
            Gradient_steps_n                    = [Gradient_steps_n; Optimisation_Values.gradient'];
            Objective_value_steps             = [Objective_value_steps; Optimisation_Values.fval];
        end

        % Geometry variable retrieval from the vector
        function [cylinder_radius, cylinder_direction_r, cylinder_direction, cylinder_centre_r, cylinder_centre] = Geometry_Retrieval(inf_cylinder_geometry_n, geometry_LB, geometry_UB, gimbal_rotation_matrix, point_cloud_centroid, Scanner_loc_cell_r, Geometry_Indices)
            % Un-normalisation
            inf_cylinder_geometry = inf_cylinder_geometry_n .* (geometry_UB - geometry_LB) + geometry_LB;
            
            % Geometry indices
            [centre_ind, azim_ind, elev_ind, radius_ind] = deal(Geometry_Indices.centre, Geometry_Indices.azim_angle, Geometry_Indices.elev_angle, Geometry_Indices.radius);

            % The spherical angles and resulting cylinder axis
            [azim_angle_r, elev_angle_r]    = deal(inf_cylinder_geometry(azim_ind), inf_cylinder_geometry(elev_ind));
            [cylinder_direction_r, ~, ~]    = Vector_Spherical_Angle_Conversion([], azim_angle_r, elev_angle_r);
    
            cylinder_direction              = (gimbal_rotation_matrix' * cylinder_direction_r')';

            % For consistency the final component is taken to always be positive
            cylinder_direction = sign(cylinder_direction(end)) * cylinder_direction;
    
            % The circle centre is converted into the cylinder centre
            circle_centre_r                 = inf_cylinder_geometry(centre_ind);
            height_r                        = 0;
            [~, ~, cylinder_centre_r, ~]    = Circle_Cylinder_Centre_Conversion(circle_centre_r, height_r, [], cylinder_direction_r, origin, Scanner_loc_cell_r);
    
            cylinder_centre_c               = (gimbal_rotation_matrix' * cylinder_centre_r')';
            cylinder_centre                 = cylinder_centre_c + point_cloud_centroid;
    
            % The cylinder radius equals the circular cross-section's radius
            cylinder_radius                 = inf_cylinder_geometry(Geometry_Indices.radius);
        end
end