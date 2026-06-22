% This function determines the fuzzy cloud based on the given QSM by calculating the distances to each cylinder and using the nearest one

function [Point_Cloud_Distributions, nearest_cylinder_ind_list] = QSM_Point_Cloud_Uncertainty(QSM, Statistical_Values, Scanner_Parameters, Scanning_Parameters, Point_Cloud_Coord)

    %% Structure inputs %%
        % QSM geometry
        cylinder_radius_list    = QSM.cylinder.radius;
        cylinder_length_list    = QSM.cylinder.length;
        cylinder_start_matrix   = QSM.cylinder.start;
        cylinder_axis_matrix    = QSM.cylinder.axis;

        % Statistical values
        Confidence_interval     = Statistical_Values.Confidence_interval;
  
        % Scanner parameters
        Empirical_Uncertainty   = Scanner_Parameters.Empirical_Uncertainty;

        % Scanning parameters
        number_scanners         = Scanning_Parameters.number_scanners;
        Scanner_loc_cell        = Scanning_Parameters.Scanner_loc_cell;

        % Point cloud
        point_cloud_matrix      = Point_Cloud_Coord.point_cloud_matrix;
        number_points_list      = Point_Cloud_Coord.number_points_list;
        
    %% Manual inputs %%
        number_chunks           = 16;               % [-] The point cloud can be split up to make the distance calculation more palatable
        Distr_Diagnostics       = false;            % [true, false]

    %% Uncertainty for each cylinder %%
        % Cylinder properties in cell arrays
        [number_cylinders, num_dim] = size(cylinder_start_matrix);
        cylinder_start_cell         = mat2cell(cylinder_start_matrix, ones(1, number_cylinders), num_dim);
        cylinder_axis_cell          = mat2cell(cylinder_axis_matrix, ones(1, number_cylinders), num_dim);
        cylinder_radius_cell        = num2cell(cylinder_radius_list);
        cylinder_length_cell        = num2cell(cylinder_length_list);

        % Distances for the point cloud to each cylinder
        number_points = sum(number_points_list);
        chunk_indices = linspace(0, number_points, number_chunks + 1);
        chunk_indices = round(chunk_indices);

        nearest_cylinder_ind_list = zeros(1, number_points);

        for c = 1 : number_chunks
            % This chunk's point cloud
            ind_start                   = chunk_indices(c) + 1;
            ind_end                     = chunk_indices(c + 1);
            chunk_point_cloud_matrix    = point_cloud_matrix(ind_start : ind_end, :);

            Point_Cylinder_Distance_fun = @(cylinder_start, cylinder_axis, cylinder_radius, cylinder_length) Point_Cylinder_Distance(cylinder_start, cylinder_axis, cylinder_radius, cylinder_length, chunk_point_cloud_matrix);
            chunk_distance_cell         = cellfun(Point_Cylinder_Distance_fun, cylinder_start_cell, cylinder_axis_cell, cylinder_radius_cell, cylinder_length_cell, 'UniformOutput', false);
    
            % Points correspond to the nearest cylinder
            chunk_distance_matrix                           = horzcat(chunk_distance_cell{:});
            [~, chunk_nearest_cylinder_ind_list]            = min(chunk_distance_matrix, [], 2);
            nearest_cylinder_ind_list(ind_start : ind_end)  = chunk_nearest_cylinder_ind_list;

            % Progress message
            if number_chunks > 1
                current_time    = datetime('now', 'format', 'dd_HHmmss');
                progress        = c / number_chunks * 100;
                fprintf('   t = %s. %g%% of the point cloud has been associated to cylinders \n', current_time, progress);
            end
        end

        % Uncertainty of each cylinder's point cloud
        [distribution_mu_cell, distribution_axes_cell, distribution_sigmae_cell, distribution_Sigma_cell] = deal(cell(1, number_points));

        for c = 1 : number_cylinders
            % Cylinder metrics
            cylinder_start  = cylinder_start_matrix(c, :);
            cylinder_axis   = cylinder_axis_matrix(c, :);

            % Points corresponding to this cylinder
            cylinder_indices_list = find(nearest_cylinder_ind_list == c);

            [cyl_distr_mu_cell, cyl_distr_axes_cell, cyl_distr_sigmae_cell, cyl_distr_Sigma_cell] = Cylinder_Uncertainty(cylinder_start, cylinder_axis, Scanner_loc_cell, Scanner_Parameters, Confidence_interval, point_cloud_matrix, number_points_list, number_scanners, cylinder_indices_list, Empirical_Uncertainty, Distr_Diagnostics);
            distribution_mu_cell(cylinder_indices_list)     = cyl_distr_mu_cell;
            distribution_axes_cell(cylinder_indices_list)   = cyl_distr_axes_cell;
            distribution_sigmae_cell(cylinder_indices_list) = cyl_distr_sigmae_cell;
            distribution_Sigma_cell(cylinder_indices_list)  = cyl_distr_Sigma_cell;
        end

        % Structure
        Point_Cloud_Distributions = struct('distribution_mu_cell', {distribution_mu_cell}, 'distribution_axes_cell', {distribution_axes_cell}, 'distribution_sigmae_cell', {distribution_sigmae_cell}, 'distribution_Sigma_cell', {distribution_Sigma_cell}, 'number_distributions', number_points, 'number_distributions_list', number_points_list);

    %% Local functions %%
        % Uncertainty for a specific cylinder
        function [distribution_mu_cell, distribution_axes_cell, distribution_sigmae_cell, distribution_Sigma_cell] = Cylinder_Uncertainty(cylinder_start, cylinder_axis, Scanner_loc_cell, Scanner_Parameters, Confidence_interval, point_cloud_matrix, number_points_list, number_scanners, cylinder_indices_list, Empirical_Uncertainty, Distr_Diagnostics)
            % This cylinder's point cloud
            cum_number_points_list      = [0, cumsum(number_points_list)];
            cylinder_point_cloud_cell   = cell(1, number_scanners);
            cyl_number_points_list      = zeros(1, number_scanners);

            for s = 1 : number_scanners
                ind_start   = cum_number_points_list(s) + 1;
                ind_end     = cum_number_points_list(s + 1);

                scanner_cyl_indices     = cylinder_indices_list(cylinder_indices_list >= ind_start & cylinder_indices_list <= ind_end);
                cylinder_point_cloud    = point_cloud_matrix(scanner_cyl_indices, :);
                
                cylinder_point_cloud_cell{s}    = cylinder_point_cloud;
                cyl_number_points_list(s)       = size(cylinder_point_cloud, 1);
            end

            Cyl_Point_Cloud_Coord   = struct('point_cloud_cell', {cylinder_point_cloud_cell}, 'number_points_list', cyl_number_points_list);

            % Uncertainty using the analytical equations
            [sigma_radial_cell, sigma_prop_cell, incidence_angle_cell, beam_range_cell, ~] = Cylindrical_Object_Uncertainty(cylinder_start, cylinder_axis, Scanner_loc_cell, Scanner_Parameters, Cyl_Point_Cloud_Coord);

            % Uncertainty derived from the empirical equation
            if Empirical_Uncertainty == true
                sigma_prop_cell = cell(1, number_scanners);

                for s = 1 : number_scanners
                    % The radial uncertainty, incidence angles and range of each point
                    sigma_radial_list       = sigma_radial_cell{s};
                    incidence_angle_list    = incidence_angle_cell{s};
                    range_list              = beam_range_cell{s};

                    % Propagation uncertainty
                    sigma_0_list        = 4e-7*range_list.^2 + 1e-5*range_list + 4e-4;
                    sigma_prop_list     = sigma_0_list + sigma_radial_list .* tan(incidence_angle_list);
                    sigma_prop_cell{s}  = sigma_prop_list;
                end
            end
    
            % The 3D distributions
            alpha                       = 1 - Confidence_interval/100;          % The confidence interval is changed to alpha
            range_bias                  = Scanner_Parameters.range_bias;
            Cylinder_Distributions      = Point_Cloud_Multivariate_Normal_Generation(alpha, Cyl_Point_Cloud_Coord, sigma_radial_cell, sigma_prop_cell, Scanner_loc_cell, range_bias, Distr_Diagnostics);

            distribution_axes_cell      = Cylinder_Distributions.distribution_axes_cell;
            distribution_sigmae_cell    = Cylinder_Distributions.distribution_sigmae_cell;
            distribution_mu_cell        = Cylinder_Distributions.distribution_mu_cell;
            distribution_Sigma_cell     = Cylinder_Distributions.distribution_Sigma_cell;

            % The covariance matrices are ensured to be symmetric as very minor asymmetry may be present
            Symmetry_fun            = @(distribution_Sigma) (distribution_Sigma + distribution_Sigma') / 2;
            distribution_Sigma_cell = cellfun(Symmetry_fun, distribution_Sigma_cell, 'UniformOutput', false);    
        end
end