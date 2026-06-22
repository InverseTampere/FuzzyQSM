% The expected Mahalanobis distance is approximated for a cylinder by taking the tangent line for each distribution
% The distance moment may be 1 or 2

function [expected_Mahal_distance, expected_Mahal_dist_list] = Expected_Mahal_Distance_Cylinder_Line_Approx(cylinder_centre, cylinder_radius, cylinder_direction, distance_moment, Point_Cloud_Distributions, Scanner_loc_cell, Plot)
    
    %% Manual inputs %%
        Point_Distance_Measure      = 'Mahalanobis';        % [Mahalanobis, Euclidean] The Euclidean or Mahalanobis distance can be used for the distance between mu and any points on the circle
        Envelope_Sampling           = true;                 % [true, false] The tangent lines of samples across part of the circle can be used to determine the weighted average E[M^m]
                                                            %               Otherwise only the nearest point (in terms of Mahal. or Euclid. distance) is used
        number_samples              = 1e3;                  % [-] If full circle sampling is used, this number of samples is taken
        Proj_Diagnostics            = false;                % [true, false] Projection of the distributions onto the cross-sectional plane

    %% Point cloud projection onto cross-section %%
        % Projection s.t. the propagation axis is always (approximately) in the -y direction
        [mu_matrix, sigmae_matrix] = Distribution_Cylinder_Cross_Section_Projection(Point_Cloud_Distributions, Scanner_loc_cell, cylinder_centre, cylinder_direction);

        % Division into dimensions for clarity
        [mu_x_list, mu_y_list]          = Column_Deal(mu_matrix);
        [sigma_x_list, sigma_y_list]    = Column_Deal(sigmae_matrix);

    %% Expected Mahalanobis distance %%
        % Samples are taken across the circle
        if Envelope_Sampling == true
            % Samples on the circle
            dtheta                  = 2*pi / number_samples;
            theta_list              = 0 : dtheta : 2*pi - dtheta;           % This avoids having duplicate samples at 0 and 2pi
            circle_sample_matrix    = cylinder_radius * [cos(theta_list)', sin(theta_list)'];

            % The expected (squared) Mahalanobis distances of each distribution
            Expected_Mahal_Distance_fun = @(mu_x, mu_y, sigma_x, sigma_y) Expected_Mahal_Distance_Circle_Line_Approx([mu_x, mu_y], [sigma_x, sigma_y], circle_sample_matrix, distance_moment, Point_Distance_Measure);
            expected_Mahal_dist_list    = arrayfun(Expected_Mahal_Distance_fun, mu_x_list, mu_y_list, sigma_x_list, sigma_y_list);

        % Only the nearest point in a Euclidean or Mahalanobis sense is used
        else
            %  The Euclidean-nearest points
            if strcmp(Point_Distance_Measure, 'Euclidean')
                % The points defining the tangent line lie on the vector towards mu
                mu_norm_list        = sqrt(sum(mu_matrix.^2, 2));                                       %#ok<*UNRCH> Unreachable code is intentional due to Point_Distance_Measure and therefore suppressed 
                line_point_matrix   = cylinder_radius * mu_matrix ./ mu_norm_list;
    
            % The Mahalanobis-nearest points
            elseif strcmp(Point_Distance_Measure, 'Mahalanobis')
                % Transformation s.t. the distribution is standard-normal
                mu_matrix_SN    = mu_matrix ./ sigmae_matrix;
                ellipse_radii   = cylinder_radius./sigmae_matrix;
                ellipse_axes    = eye(2);
                ellipse_centre  = zeros(1, 2);
    
                % Projection onto the ellipse leads to the nearest points in the transformed coordinate frame which define the tangent lines
                Projection_fun          = @(radius_x, radius_y, x, y) Point_to_Ellipse_Projection(ellipse_centre, [radius_x, radius_y], ellipse_axes, [x, y], Proj_Diagnostics, Proj_Diagnostics);
                line_point_SN_cell      = arrayfun(Projection_fun, ellipse_radii(:, 1), ellipse_radii(:, 2), mu_matrix_SN(:, 1), mu_matrix_SN(:, 2), 'UniformOutput', false);
                line_point_SN_matrix    = vertcat(line_point_SN_cell{:});
                line_point_matrix       = line_point_SN_matrix .* sigmae_matrix;
            end

            % The expected (squared) Mahalanobis distances of each distribution
            [p_x_list, p_y_list]        = Column_Deal(line_point_matrix);       % For clarity the line points are divided into their dimensions

            Expected_Mahal_Distance_fun = @(mu_x, mu_y, sigma_x, sigma_y, p_x, p_y) Expected_Mahal_Distance_Circle_Line_Approx([mu_x, mu_y], [sigma_x, sigma_y], [p_x, p_y], distance_moment, Point_Distance_Measure);
            expected_Mahal_dist_list    = arrayfun(Expected_Mahal_Distance_fun, mu_x_list, mu_y_list, sigma_x_list, sigma_y_list, p_x_list, p_y_list);
        end

        % The average over all distributions
        expected_Mahal_distance = mean(expected_Mahal_dist_list);

    %% Plot %%
        if Plot == true
            % Distribution properties
            total_num_distributions     = Point_Cloud_Distributions.number_distributions;
            distribution_mu_cell        = Point_Cloud_Distributions.distribution_mu_cell;
            distribution_axes_cell      = Point_Cloud_Distributions.distribution_axes_cell;
            distribution_sigmae_cell    = Point_Cloud_Distributions.distribution_sigmae_cell;

            % The number of coordinates used for each distribution and the cylinder
            number_coord = 1e2;     

            % Estimate of the length for the plot
            mu_0_matrix         = vertcat(distribution_mu_cell{:}) - cylinder_centre;
            mu_0_norm_list      = sqrt(sum(mu_0_matrix.^2, 2));
            cylinder_length     = 2*max(mu_0_norm_list);

            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            hold on
            grid on
            
            % Cylinder
            [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, cylinder_length, cylinder_centre, cylinder_direction, number_coord);
            surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Cylinder');

            % Distributions and their distances           
            for d = 1 : total_num_distributions
                % Distribution
                proj_sigmae = sigmae_matrix(d, :);          % Used to scale the Mahalanobis distance

                mu          = distribution_mu_cell{d};
                sigmae      = distribution_sigmae_cell{d};
                distr_axes  = distribution_axes_cell{d};

                [distr_coord_matrix, number_coord_distr] = Ellipsoid_Coordinate_Generator(mu, sigmae, distr_axes, number_coord);

                x_ellipsoid = reshape(distr_coord_matrix(:, 1), sqrt(number_coord_distr) * [1, 1]);
                y_ellipsoid = reshape(distr_coord_matrix(:, 2), sqrt(number_coord_distr) * [1, 1]);
                z_ellipsoid = reshape(distr_coord_matrix(:, 3), sqrt(number_coord_distr) * [1, 1]);
            
                surf_distr  = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.10, 'DisplayName', '1 \sigma');
                sc_mu       = scatter3(mu(1), mu(2), mu(3), 'filled', 'MarkerFaceColor', 'r', 'DisplayName', '\mu');
                
                % Mahalanobis distance, adjusted to be in metres
                Mahal_dist  = expected_Mahal_dist_list(d)^(1/distance_moment);

                [proj_mu, ~, ~] = Point_to_Vector_Projection(mu, cylinder_direction, cylinder_centre);

                vector          = mu - proj_mu;

                vector_sign     = sign(cylinder_radius - norm(vector));     % Positive if it points outward
                vector          = sqrt(prod(proj_sigmae)) * Mahal_dist * vector_sign * vector / norm(vector);
                
                if distance_moment == 1
                    vec_name = sprintf('E[M] * sqrt(%s %s)', '\sigma_x', '\sigma_y');
                elseif distance_moment == 2
                    vec_name = sprintf('sqrt(E[M^2] * %s %s)', '\sigma_x', '\sigma_y');
                end
                
                pl_vec = plot3(mu(1) + [0, vector(1)], mu(2) + [0, vector(2)], mu(3) + [0, vector(3)], 'color', 'k', 'lineWidth', 1, 'DisplayName', vec_name);
                
                if d > 1
                    surf_distr.HandleVisibility = 'Off';
                    sc_mu.HandleVisibility      = 'Off';
                    pl_vec.HandleVisibility     = 'Off';
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
            
            fprintf('E[M^%g] = %.3g for the cylinder using the line approximation. \n', distance_moment, expected_Mahal_distance);
            disp('The script will finish and figure will close upon a key-press.');
            pause();
            
            close(1);
        end

    %% Local functions %%
        % The weighted-average expected Mahalanobis distance is computed over the tangent lines located on the given points
        function expected_Mahal_distance = Expected_Mahal_Distance_Circle_Line_Approx(distr_mu, distr_sigmae, line_point_matrix, distance_moment, Weight_Distance_Measure)
            % Individual dimensions for clarity
            [mu_x, mu_y]                    = Column_Deal(distr_mu);
            [sigma_x, sigma_y]              = Column_Deal(distr_sigmae);
            [x_circle_list, y_circle_list]  = Column_Deal(line_point_matrix);
            
            % Expected (squared) Mahalanobis distance to tangent lines placed at each point
            tau_list = abs(x_circle_list.*(x_circle_list - mu_x) + y_circle_list.*(y_circle_list - mu_y)) ./ sqrt(sigma_x^2*x_circle_list.^2 + sigma_y^2*y_circle_list.^2);
    
            if distance_moment == 1
                point_expected_Mahal_dist_list = tau_list.*erf(tau_list/sqrt(2)) + sqrt(2/pi)*exp(-1/2*tau_list.^2);
            elseif distance_moment == 2
                point_expected_Mahal_dist_list = tau_list.^2 + 1;
            else
                error('The distance moment can only be 1 or 2.');
            end
    
            % Each point's weight follows from the distance between mu and the point
            if strcmp(Weight_Distance_Measure, 'Euclidean')
                distance_list = sqrt((x_circle_list - mu_x).^2 + (y_circle_list - mu_y).^2);
            elseif strcmp(Weight_Distance_Measure, 'Mahalanobis')
                distance_list = sqrt((x_circle_list - mu_x).^2/sigma_x^2 + (y_circle_list - mu_y).^2/sigma_y^2);
            end
    
            weight_list = 1./distance_list.^(distance_moment + 1);      % The inverse is used as a greater distance is worse. 
                                                                        % The distance moment is increased to weigh this effect more then the expected Mahalanobis distance
            weight_list = weight_list / sum(weight_list);               % s.t. the sum is 1
    
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
end