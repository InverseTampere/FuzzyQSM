% This script takes an initial cylinder and its point cloud and fits a fuzzy circle cross-section to it
% The cylidner axis thus remains the same

function [Fuzzy_Inf_Cylinder, Point_Cloud_Distributions] = Fuzzy_QSM_Circle_Fitting(Init_Cylinder, Point_Cloud_Coord, Fitting_Parameters, Scanner_Parameters, Scanning_Parameters)

    %% Structure inputs %%
        % Initial cylinder
        init_cylinder_start     = double(Init_Cylinder.start);          % Converted to double for fmincon
        init_cylinder_axis      = double(Init_Cylinder.axis);
        init_cylinder_radius    = double(Init_Cylinder.radius);
        init_cylinder_length    = double(Init_Cylinder.length);

        % Fitting Parameters
        bounds_margin           = Fitting_Parameters.bounds_margin;

        % Scanner parameters
        range_bias              = Scanner_Parameters.range_bias;

        % Scanning parameters
        Scanner_loc_cell        = Scanning_Parameters.Scanner_loc_cell;
        number_scanners         = Scanning_Parameters.number_scanners;

    %% Manual inputs %%
        % If desired, intermediate results can be checked
        Distr_Diagnostics       = false;        % [true, false]
        Proj_Diagnostics        = false;        % [true, false]

        % Optimiser inputs
        max_iterations          = 1e2;          % [-] Maximum number of iterations by the optimiser. If it is reached, a warning is displayed
        function_tolerance      = 1e-4;         % [-] Minimum size for the step in the objective function. Note that it is normalised by the initial value
        step_tolerance          = 1e-4;         % [-] Minimum size for the step in the design parameters. Note that they are normalised between [0, 1]
        FD_step_size            = 1e-4;         % [-] Step size taken to evaluate the finite difference. Note that sqrt(eps) is the default
       
        % Objective function inputs
        number_samples          = 1e3;          % [-] Number of samples around the circle
        distance_moment         = 2;            % [-] The Mahalanobis distance is taken to this power

        % Displayed results
        Print                   = false;         % [true, false]
        Plot                    = false;         % [true, false]

    %% Fuzzy cloud %%
        % Uncertainty based on the initial cylinder geometry
        init_cylinder_centre                            = init_cylinder_start + init_cylinder_length/2 * init_cylinder_axis;
        [sigma_radial_cell, sigma_prop_cell, ~, ~, ~]   = Cylindrical_Object_Uncertainty(init_cylinder_centre, init_cylinder_axis, Scanner_loc_cell, Scanner_Parameters, Point_Cloud_Coord);

        % The 3D distributions
        alpha                       = 0.99;          % The alpha value is irrelevant
        Point_Cloud_Distributions   = Point_Cloud_Multivariate_Normal_Generation(alpha, Point_Cloud_Coord, sigma_radial_cell, sigma_prop_cell, Scanner_loc_cell, range_bias, Distr_Diagnostics);

        % Projection onto the cross-section
        [~, cyl_vector_basis, ~]    = Vector_Based_Rotation(init_cylinder_centre, init_cylinder_axis, init_cylinder_centre);
        Projected_Distributions     = Multivariate_Normal_Plane_Projection(cyl_vector_basis, init_cylinder_centre, Point_Cloud_Distributions, Proj_Diagnostics);

   %% Solution bounds %%
        % The geometry parameters' indices (as the optimiser's inputs must be vector format)
        centre_ind                  = 1 : 2;
        radius_ind                  = max(centre_ind) + 1;
        num_optim_geo_parameters    = radius_ind;

        Geometry_Indices            = struct('centre', centre_ind, 'radius', radius_ind);

        % The set of initial centered geometry parameters
        num_dim                                         = length(init_cylinder_axis);
        init_circle_centre                              = zeros(1, num_dim - 1);                                    % Projection is performed centered on the initial cylinder start
        circle_geometry_init                            = zeros(1, num_optim_geo_parameters);
        circle_geometry_init([centre_ind, radius_ind])  = [init_circle_centre, init_cylinder_radius];       

        % Reasonable bounds for the optimisation process are based on the least-squares geometry s.t. the normalised variables can be directly compared
        circle_centre_LB = init_circle_centre - bounds_margin * init_cylinder_radius;
        circle_centre_UB = init_circle_centre + bounds_margin * init_cylinder_radius;
        
        radius_LB = init_cylinder_radius / (1 + bounds_margin);
        radius_UB = init_cylinder_radius * (1 + bounds_margin);

        % The combined sets of bounds
        geometry_LB                             = zeros(1, num_optim_geo_parameters);
        geometry_LB([centre_ind, radius_ind])   = [circle_centre_LB, radius_LB];

        geometry_UB                             = zeros(1, num_optim_geo_parameters);
        geometry_UB([centre_ind, radius_ind])   = [circle_centre_UB, radius_UB];

        % As a result, the optimisation bounds are between 0 and 1
        Optim_LB = zeros(1, num_optim_geo_parameters);
        Optim_UB = ones(1, num_optim_geo_parameters);
        
    %% Circle %%
        %--% Optimisation %--%
        % Normalised geometry and the initial objective value for normalisation
        circle_geometry_init_n      = (circle_geometry_init - geometry_LB) ./ (geometry_UB - geometry_LB);
        circle_objective_value_init = Circle_Objective_Function(circle_geometry_init_n, Geometry_Indices, Projected_Distributions, geometry_LB, geometry_UB, 1, number_samples, distance_moment);

        Circle_Objective_Function_fun = @(circle_geometry_n) Circle_Objective_Function(circle_geometry_n, Geometry_Indices, Projected_Distributions, geometry_LB, geometry_UB, circle_objective_value_init, number_samples, distance_moment);

        % Certain information is saved during optimisation
        [circle_geometry_steps_n, Gradient_steps_n, Objective_value_steps] = deal([]);

        % Optimisation with the interior-point algorithm
        Options = optimoptions('fmincon', 'Display', 'off', 'OutputFcn', @Optim_Path, 'Algorithm', 'interior-point', 'SubProblemAlgorithm', 'cg', 'EnableFeasibilityMode', true, 'MaxIterations', max_iterations, 'FunctionTolerance', function_tolerance, 'StepTolerance', step_tolerance, 'FiniteDifferenceStepSize', FD_step_size);

        [circle_geometry_n, Objective_value, ~, Output, ~, ~, Optimum_Hessian_matrix_n] = fmincon(Circle_Objective_Function_fun, circle_geometry_init_n, [], [], [], [], Optim_LB, Optim_UB, [], Options);

        %--% Checks %--%
        % Check if the optimum is worse than the initial geometry, in which case sqp is used as the optimisation algorithm
        if Objective_value > 1
            % The optimisation data is cleared
            [circle_geometry_steps_n, Gradient_steps_n, Objective_value_steps] = deal([]);

            % Updated optimisation method
            Options.Algorithm = 'sqp';
            
            % New optimisation run
            [circle_geometry_n, Objective_value, ~, Output, ~, ~, Optimum_Hessian_matrix_n] = fmincon(Circle_Objective_Function_fun, circle_geometry_init_n, [], [], [], [], Optim_LB, Optim_UB, [], Options);    % Optimisation with sqp
        end
        
        % Check if the maximum number of iterations was reached
        number_optimiser_steps = Output.iterations;

        if number_optimiser_steps == max_iterations
            warning('optimiser:max_iter', 'The maximum number of iterations for the optimiser was reached. Check whether this should be increased or whether the objective function is poor.')
        end

        % Check if any of the bounds are reached, which may indicate that the optimum lies outside the optimisation space
        [centre_n, radius_n] = deal(circle_geometry_n(centre_ind), circle_geometry_n(radius_ind));

        if max(centre_n) > 1 - 2*step_tolerance || min(centre_n) < 2*step_tolerance                                                                             % Note that the finite step size means the bounds are unlikely to be hit exactly
            warning('bounds:centre', 'The bounds for optimisation have restricted the optimal centre location. Try increasing the bounds margin.');
        end
        if radius_n > 1 - 2*step_tolerance || radius_n < 2*step_tolerance
            warning('bounds:radius', 'The bounds for optimisation have restricted the optimal radius. Try increasing the bounds margin.');
        end

        % Check if the objective function is convex at the solution (i.e. if all eigenvalues of the Hessian are positive)
        % Note that as the optimising algorithm is quasi-Newton, convexity is not guaranteed
        Hessian_n_eigenvalues       = eig(Optimum_Hessian_matrix_n);
        Hessian_n_eigenvalues_sign  = sign(Hessian_n_eigenvalues);

        if sum(Hessian_n_eigenvalues_sign) < num_optim_geo_parameters
            warning('optimiser:convexity', 'The objective function is not convex at the found solution. The lowest Hessian eigenvalue is %.3g \n', min(Hessian_n_eigenvalues));
        end

        % Circle geometry parameters
        circle_geometry = circle_geometry_n .* (geometry_UB - geometry_LB) + geometry_LB;
        circle_centre   = circle_geometry(centre_ind);
        circle_radius   = circle_geometry(radius_ind);

        % 3D cylinder centre
        cylinder_centre_t   = (cyl_vector_basis' * [circle_centre, 0]')';
        cylinder_centre     = cylinder_centre_t + init_cylinder_centre;

        % Structure containing the infinite cylinder geometry
        Fuzzy_Inf_Cylinder = struct('radius', circle_radius, 'axis', init_cylinder_axis, 'centre', cylinder_centre);

    %% Fitting information %%
        % Printed statement
        if Print == true
            % The objective value is converted into the expected Mahalanobis distance
            Expected_Mahal_distance = Objective_value * circle_objective_value_init;

            fprintf('The infinite cylinder has been fitted with O =  %.4g, E[M^2] = %.4g \n', Objective_value, Expected_Mahal_distance);
            fprintf('   Radius:     %.3g m \n',             circle_radius);
            fprintf('   Centre:     [%.3g, %.3g, %.3g] \n', cylinder_centre);
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
            [init_cylinder_coord_x, init_cylinder_coord_y, init_cylinder_coord_z, ~] = Cylinder_Surface_Generator(init_cylinder_radius, init_cylinder_length, init_cylinder_centre, init_cylinder_axis, number_coord);
            surf(init_cylinder_coord_x, init_cylinder_coord_y, init_cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'm', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Initial cylinder');
            
            % The optimised cylinder surface
            [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(circle_radius, init_cylinder_length, cylinder_centre, init_cylinder_axis, number_coord);
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
            Geometry_steps_n_delta_matrix       = diff(circle_geometry_steps_n, 1, 1);
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
        % Circle fitting objective function
        function Objective_value = Circle_Objective_Function(circle_geometry_n, Geometry_Indices, Projected_Distributions, geometry_LB, geometry_UB, Objective_value_init, number_samples, distance_moment)

            %% Structure inputs %%
                % Geometry indices
                centre_indices          = Geometry_Indices.centre;
                radius_index            = Geometry_Indices.radius;

                % Projected distributions
                proj_mu_cell            = Projected_Distributions.Plane.Projection.mu;
                proj_sigmae_cell        = Projected_Distributions.Plane.Projection.sigmae;
                proj_axes_cell          = Projected_Distributions.Plane.Projection.distr_axes;

            %% Manual inputs %%
                % It is recommend that outputs are set to false        
                Objective_Diagnostics   = false;        % [true, false] Shows the objective values for each distribution

            %% Geometry parameters (optimisation variables) %%
                % Un-normalising the geometry vector
                circle_geometry = circle_geometry_n .* (geometry_UB - geometry_LB) + geometry_LB;
        
                % Circular cross section
                circle_centre   = circle_geometry(centre_indices);               
                circle_radius   = circle_geometry(radius_index);    

            %% Weighted expected Mahalanobis distance %%
                % Samples on the circle
                dtheta                  = 2*pi / number_samples;
                theta_list              = 0 : dtheta : 2*pi - dtheta;           % This avoids having duplicate samples at 0 and 2pi
                circle_sample_matrix    = circle_centre + circle_radius * [cos(theta_list)', sin(theta_list)'];

                % Weighted expected Mahalanobis distance for each distribution
                Expected_Mahal_Distance_fun = @(distr_mu, distr_sigmae, distr_axes) Expected_Mahal_Distance_Circle_Line_Approx(distr_mu, distr_sigmae, distr_axes, circle_sample_matrix, distance_moment);
                expected_Mahal_dist_list    = cellfun(Expected_Mahal_Distance_fun, proj_mu_cell, proj_sigmae_cell, proj_axes_cell);

                % The average over all distributions divided by the initial value is taken as the objective value
                distr_objective_list    = expected_Mahal_dist_list / Objective_value_init;
                Objective_value         = mean(distr_objective_list);

            %% Diagnostics %%
                % For diagnostic purposes the objective value and geometry are printed
                if Objective_Diagnostics == true
                    % Printed messages
                    disp('---');
                    fprintf('O = %.3g \n', Objective_value);
                    fprintf('c = [%.3g, %.3g] m \n', circle_centre);
                    fprintf('r = %.3g m \n', circle_radius);       
        
                    % Objective colour map
                    number_colours  = 1e3;
                    objective_cmap = cbrewer('seq', 'Blues', number_colours);
                    objective_cmap = max(objective_cmap, 0);
                    objective_cmap = min(objective_cmap, 1);
        
                    min_objective       = min(distr_objective_list);
                    max_objective       = max([distr_objective_list; min_objective + 1e-16]);                             % In case the values are all identical
                    objective_list_n    = (distr_objective_list - min_objective) / (max_objective - min_objective);
                    cmap_ind            = round((number_colours - 1) * objective_list_n) + 1;
        
                    % Number of coordinates per circle/distribution
                    number_coordinates = 1e2;
        
                    % Plot
                    figure(1)
                    % Set the size and white background color
                    set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                    set(gcf, 'color', [1, 1, 1])    
        
                    % First axes for the cylinder and distribution
                    ax1 = axes;
                    hold on
                    grid on
        
                    % Circle
                    circle_coord_matrix = Circle_Coordinates(circle_centre, circle_radius, [], number_coordinates);
                    pl_circle_name      = 'Circle';
                    pl_circle           = plot(ax1, circle_coord_matrix(:, 1), circle_coord_matrix(:, 2), 'LineWidth', 2, 'Color', 'b', 'DisplayName', pl_circle_name);

                    % Distributions
                    number_distributions = length(proj_mu_cell);

                    for i = 1 : number_distributions
                        % Distribution properties
                        distr_mu        = proj_mu_cell{i};
                        distr_axes      = proj_axes_cell{i};
                        distr_sigmae    = proj_sigmae_cell{i};
        
                        % Coordinates
                        distr_coord_matrix = Ellipse_Coordinate_Generator(distr_mu, distr_axes, distr_sigmae, number_coordinates);
            
                        % Distribution at 1 sigma
                        pl_distr_name   = sprintf('Distribution, 1%s', '\sigma');
                        pl_distr        = plot(distr_coord_matrix(:, 1), distr_coord_matrix(:, 2), 'LineWidth', 1, 'Color', 'r', 'DisplayName', pl_distr_name);
                    end
        
                    % Second axes for the points and their objective values
                    ax2 = axes;
                    hold on
                    grid on
        
                    proj_mu_matrix = vertcat(proj_mu_cell{:});
                    sc_points_name = 'Point cloud';
                    sc_points      = scatter(ax2, proj_mu_matrix(:, 1), proj_mu_matrix(:, 2), 200, objective_cmap(cmap_ind, :), 'Marker', '.', 'DisplayName', sc_points_name);
        
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
                
                    ax1.DataAspectRatio = [1, 1, 1];
                    ax2.DataAspectRatio = [1, 1, 1];
                        
                    % Legend
                    legend(ax1, [pl_circle, pl_distr, sc_points], {pl_circle_name, pl_distr_name, sc_points_name}, 'location', 'northoutside');
            
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
        
                    % Pause message
                    disp('The objective values have been determined. The figure will close and script end upon a key-press.');
                    pause();
            
                    close(1);
                end
        end
    
        % The weighted-average expected Mahalanobis distance is computed over the tangent lines located on the given points
        function expected_Mahal_distance = Expected_Mahal_Distance_Circle_Line_Approx(distr_mu, distr_sigmae, distr_axes, circle_sample_matrix, distance_moment)
            % Transformation to the distribution's coordinate frame
            circle_sample_matrix_t = circle_sample_matrix - distr_mu;
            circle_sample_matrix_r = (distr_axes * circle_sample_matrix_t')';

            % Individual dimensions for clarity
            [mu_x, mu_y]                    = deal(0, 0);                               % As it is now centered on the distribution
            [sigma_x, sigma_y]              = Column_Deal(distr_sigmae);
            [x_circle_list, y_circle_list]  = Column_Deal(circle_sample_matrix_r);
            
            % Expected (squared) Mahalanobis distance to tangent lines placed at each point
            tau_list = abs(x_circle_list.*(x_circle_list - mu_x) + y_circle_list.*(y_circle_list - mu_y)) ./ sqrt(sigma_x^2*x_circle_list.^2 + sigma_y^2*y_circle_list.^2);
    
            if distance_moment == 1
                point_expected_Mahal_dist_list = tau_list.*erf(tau_list/sqrt(2)) + sqrt(2/pi)*exp(-1/2*tau_list.^2);
            elseif distance_moment == 2
                point_expected_Mahal_dist_list = tau_list.^2 + 1;
            else
                error('The distance moment can only be 1 or 2.');
            end
    
            % Each point's weight follows from the Mahalanobis distance between mu and the point
            distance_list   = sqrt((x_circle_list - mu_x).^2/sigma_x^2 + (y_circle_list - mu_y).^2/sigma_y^2);
            weight_list     = 1./distance_list.^(distance_moment + 1);      % The inverse is used as a greater distance is worse. 
                                                                            % The distance moment is increased to weigh this effect more then the expected Mahalanobis distance
            weight_list     = weight_list / sum(weight_list);               % s.t. the sum is 1
    
            % If the distance is 0 for one point it means it lies exactly at a circle point. Its weight is set to 1, with the others 0
            numerical_margin    = 1e-6;
            zero_dist_bool      = distance_list < numerical_margin;
    
            if sum(zero_dist_bool) == 1
                weight_list(zero_dist_bool)     = 1;
                weight_list(~zero_dist_bool)    = 0;
            end
    
            % The weighted expected distance is then
            expected_Mahal_distance = sum(weight_list .* point_expected_Mahal_dist_list);
        end

        % Optimisation path function
        function stop = Optim_Path(circle_geometry_n, Optimisation_Values, ~)
            stop = false;

            % Concatenate objective value, geometry and various other terms computed during each iteration
            circle_geometry_steps_n = [circle_geometry_steps_n; circle_geometry_n];
            Gradient_steps_n        = [Gradient_steps_n; Optimisation_Values.gradient'];
            Objective_value_steps   = [Objective_value_steps; Optimisation_Values.fval];
        end
end