% Multivariate normal distributions in the form of 3D ellipsoids are generated from the point cloud and its uncertainty

function Point_Cloud_Distributions = Point_Cloud_Multivariate_Normal_Generation(alpha, Point_Cloud_Coord, sigma_radial_cell, sigma_prop_cell, Scanner_loc_cell, range_bias, Plot)

    %% Inputs %%
        point_cloud_cell    = Point_Cloud_Coord.point_cloud_cell;
        number_points_list  = Point_Cloud_Coord.number_points_list;

    %% Expected values %%
        % The expected values are the point cloud locations, adjusted for range bias
        point_cloud_cell_rb = Range_Bias_Compensation_3D(point_cloud_cell, range_bias, Scanner_loc_cell); 

        mu_matrix                       = vertcat(point_cloud_cell_rb{:});
        [num_distributions, num_dim]    = size(mu_matrix);
        distribution_mu_cell            = mat2cell(mu_matrix, ones(1, num_distributions), num_dim);
        
    %% Uncertainty and coverage probability %%
        % The radial uncertainty acts on two of the distribution's axes, the range uncertainty on one
        sigma_radial_list           = vertcat(sigma_radial_cell{:});
        sigma_prop_list             = vertcat(sigma_prop_cell{:});
        distribution_sigmae_matrix  = [sigma_radial_list, sigma_radial_list, sigma_prop_list];
        distribution_sigmae_cell    = mat2cell(distribution_sigmae_matrix, ones(1, num_distributions), num_dim);

        % The radii are these uncertainties, taking coverage probability into account
        m_STD                       = sqrt(chi2inv(1 - alpha, num_dim));          % The number of standard deviations in each direction to achieve alpha
        
        distribution_radii_matrix   = m_STD * distribution_sigmae_matrix;
        distribution_radii_cell     = mat2cell(distribution_radii_matrix, ones(1, num_distributions), num_dim);

    %% Distribution orientation %%
        % The axes are dependent on the scanner location
        num_scanners            = length(Scanner_loc_cell);
        distr_scanner_loc_cell  = cell(1, num_scanners);
        distribution_axes_cell  = cell(1, num_scanners);

        for s = 1 : num_scanners
            % This scanner's point location and point cloud
            scanner_loc     = Scanner_loc_cell{s};
            points_matrix   = point_cloud_cell_rb{s};
            num_points      = number_points_list(s);
            
            % Vectors from scanner to points, i.e. the propagation axis
            prop_axis_matrix        = points_matrix - scanner_loc;
            prop_axis_matrix        = prop_axis_matrix ./ sqrt(sum(prop_axis_matrix.^2, 2));
                                
            % The radial axes are orthogonal to their respective range axis
            radial_axis_1_matrix    = [prop_axis_matrix(:, 2), -prop_axis_matrix(:, 1), zeros(num_points, 1)];
            radial_axis_2_matrix    = cross(prop_axis_matrix, radial_axis_1_matrix, 2);
            
            % The combined axes
            prop_axis_matrix        = reshape(prop_axis_matrix', [1, num_dim, num_points]);
            radial_axis_1_matrix    = reshape(radial_axis_1_matrix', [1, num_dim, num_points]);
            radial_axis_2_matrix    = reshape(radial_axis_2_matrix', [1, num_dim, num_points]);
            
            axis_matrix             = [radial_axis_1_matrix; radial_axis_2_matrix; prop_axis_matrix];
            
            axis_cell                   = mat2cell(axis_matrix, num_dim, num_dim, ones(1, num_points));
            distribution_axes_cell{s}   = squeeze(axis_cell);

            % The scanner locations are repeated
            distr_scanner_loc_cell{s}  = repmat({scanner_loc}, num_points, 1);
        end
        
        % Combined axes
        distribution_axes_cell          = vertcat(distribution_axes_cell{:});
        distribution_scanner_loc_cell   = vertcat(distr_scanner_loc_cell{:});
        
    %% Covariance matrices %%
        % The covariance matrix follows from rotating the diagonal uncertainty matrix by the inverse of its axes
        Covariance_fun          = @(sigmae, axes) axes \ (sigmae.^2 .* eye(num_dim)) * axes;
        distribution_Sigma_cell = cellfun(Covariance_fun, distribution_sigmae_cell, distribution_axes_cell, 'UniformOutput', false);
        
        % The covariance matrices are ensured to be symmetric as very minor asymmetry may be present
        Symmetry_fun            = @(distribution_Sigma) (distribution_Sigma + distribution_Sigma') / 2;
        distribution_Sigma_cell = cellfun(Symmetry_fun, distribution_Sigma_cell, 'UniformOutput', false);    

    %% Distributions structure %%
        % The uncertainty and distributions are saved in a structure
        Point_Cloud_Distributions = struct('number_distributions_list', number_points_list, 'number_scanners', length(Scanner_loc_cell), 'number_distributions', num_distributions, 'sigma_radial_cell', {num2cell(sigma_radial_list)}, 'sigma_prop_cell', {num2cell(sigma_prop_list)}, ...
                                           'distribution_axes_cell', {distribution_axes_cell}, 'distribution_sigmae_cell', {distribution_sigmae_cell}, 'alpha', alpha, 'distribution_radii_cell', {distribution_radii_cell}, 'distribution_mu_cell', {distribution_mu_cell}, 'distribution_Sigma_cell', {distribution_Sigma_cell}, 'distribution_scanner_loc_cell', {distribution_scanner_loc_cell});
        
    %% Plot %%
        % Plot showing all the lasers, points and their multivariate distributions (ellipsoids)
        num_coords = 1e2;
        
        if Plot == true
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            grid on
            hold on
            
            % Ellipsoids
            for e = 1: num_distributions
                % Ellipsoid properties
                ellipsoid_centre    = distribution_mu_cell{e};
                ellipsoid_radii     = distribution_radii_matrix(e, :);
                ellipsoid_axes      = distribution_axes_cell{e};
                
                % Surface
                [ellipsoid_coord_matrix, num_coords] = Ellipsoid_Coordinate_Generator(ellipsoid_centre, ellipsoid_radii, ellipsoid_axes, num_coords);

                x_ellipsoid = reshape(ellipsoid_coord_matrix(:, 1), sqrt(num_coords) * [1, 1]);
                y_ellipsoid = reshape(ellipsoid_coord_matrix(:, 2), sqrt(num_coords) * [1, 1]);
                z_ellipsoid = reshape(ellipsoid_coord_matrix(:, 3), sqrt(num_coords) * [1, 1]);

                surf_el = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'DisplayName', sprintf('Multivariate normal, \\alpha = %.3g', alpha));

                % Axes
                for d = 1 : num_dim
                    ellipsoid_axis = ellipsoid_radii(d) * ellipsoid_axes(d, :);   
                    ax_el = plot3(ellipsoid_centre(1) + [0, ellipsoid_axis(1)], ellipsoid_centre(2) + [0, ellipsoid_axis(2)], ellipsoid_centre(3) + [0, ellipsoid_axis(3)], 'LineWidth', 0.5, 'color', 'k', 'DisplayName', 'Ellipsoid axes');

                    if e > 1 || d > 1
                        ax_el.HandleVisibility = 'Off';
                    end                    
                end
                
                % Laser beam
                cum_points_list = cumsum(number_points_list);
                diff_list       = cum_points_list - e;
                scanner_ind     = diff_list == min(diff_list(diff_list >= 0));
                scanner_loc     = Scanner_loc_cell{scanner_ind};
                
                pl_las = plot3([scanner_loc(1), ellipsoid_centre(1)], [scanner_loc(2), ellipsoid_centre(2)], [scanner_loc(3), ellipsoid_centre(3)], 'LineWidth', 0.5, 'color', 'r', 'DisplayName', 'Laser');
                
                if e > 1
                    surf_el.HandleVisibility    = 'Off';
                    pl_las.HandleVisibility     = 'Off';
                end
            end
            
            % Axis labels
            xlabel('x [m]');
            ylabel('y [m]');
            zlabel('z [m]');
            
            % Axis bounds
            axis equal   
            
            largest_radius  = max(distribution_radii_matrix, [], 'all');
            
            ellipsoid_LB    = min(mu_matrix, [], 1) - largest_radius; 
            ellipsoid_UB    = max(mu_matrix, [], 1) + largest_radius;
            
            xlim([ellipsoid_LB(1), ellipsoid_UB(1)]);
            ylim([ellipsoid_LB(2), ellipsoid_UB(2)]);
            zlim([ellipsoid_LB(3), ellipsoid_UB(3)]);
                        
            view(45, 45);

            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    
            
            % Pause
            disp('The script will continue and the figure will close when a key is pressed');
            pause();
            close(1);
        end
    
end