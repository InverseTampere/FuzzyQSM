% This script takes a noise-free Helios scan and adds Gaussian noise to it to simulate an actual laser scan
% Note that it also removes the points which exceed the maximum incidence angle

% Inputs:
%   Point_Cloud_Data_n:
%       point_cloud_cell                : [m] Cell array containing the noise-free point cloud of each scanner {[n1 x 3], [n2 x 3], ...}
%       number_points_list              : [-] Array of the number of points of each scanner [n1, n2, ...]
%
%   Scanner_Parameters:
%       beam divergence                 : [rad] 2 sigma (1/e2) divergence half-angle of the beam
%       beam_exit_diameter              : [m] Diameter of the beam when exiting the scanner
%       sigma_range_0                   : [m] 1 sigma range uncertainty at zero incidence angle
%       range_bias                      : [m] Range bias of the scanner
%       max_incidence_angle             : [rad] Maximum allowed incidence angle
%
%   Scanning_Parameters:
%       number_scanners                 : [-] Number of given scanner locations
%       Scanner_loc_cell                : [m] Scanner locations {[x, y, z], [x, y, z], ...}
%
%   Triangular_Mesh:
%       Triangles.normal_vector         : [-] Vectors normal to each triangle plane, [t x 3] where t is the number of triangles
%
%   triangle_index_cell                 : [-] Cell array containing the triangle index corresponding to each point {[n1 x 1], [n2 x 1], ...}

% Outputs:
%   Point_Cloud_Data_n:
%       point_cloud_cell                : [m] Cell array containing the noisy point cloud of each scanner {[n1 x 3], [n2 x 3], ...}
%       number_points_list              : [-] Array of the number of points of each scanner. Note that it is not necessarily the same as the input. [n1, n2, ...]

function Point_Cloud_Data_n = Helios_Scan_Noise_Addition(Point_Cloud_Data, Scanner_Parameters, Scanning_Parameters, Triangular_Mesh, triangle_index_cell, Plot)

    %% Structure inputs %%
        % Point cloud data
        point_cloud_cell                = Point_Cloud_Data.point_cloud_cell;
        number_points_list              = Point_Cloud_Data.number_points_list;

        % Scanner parameters
        beam_divergence                 = Scanner_Parameters.beam_divergence;
        beam_exit_diameter              = Scanner_Parameters.beam_exit_diameter;          
        sigma_range_0                   = Scanner_Parameters.sigma_range_0;       
        range_bias                      = Scanner_Parameters.range_bias;
        max_incidence_angle             = Scanner_Parameters.max_incidence_angle;

        % Scanning parameters
        Scanner_loc_cell                = Scanning_Parameters.Scanner_loc_cell;
        number_scanners                 = Scanning_Parameters.number_scanners;

        % Triangular mesh
        triangle_normal_vector_matrix   = Triangular_Mesh.Triangles.normal_vector;
    
    %% The noisy point cloud %%
        point_cloud_cell_n = cell(1, number_scanners);

        for s = 1 : number_scanners
            % This scanner's data
            scanner_location    = Scanner_loc_cell{s};
            point_cloud_matrix  = point_cloud_cell{s};
            triangle_index_list = triangle_index_cell{s};
            number_points       = number_points_list(s);

            % The vectors to the points
            vector_matrix       = point_cloud_matrix - scanner_location;
            point_range_list    = sqrt(sum(vector_matrix.^2, 2));
            vector_matrix       = vector_matrix ./ point_range_list;

            % The normal vectors
            normal_vector_matrix = triangle_normal_vector_matrix(triangle_index_list, :);

            % Resulting incidence angles
            incidence_angle_list = acos(abs(dot(vector_matrix, normal_vector_matrix, 2)));

            % Range uncertainty 
            beam_width_list     = beam_exit_diameter + 2*point_range_list * tan(beam_divergence);
            sigma_radial_list   = beam_width_list / 4;                                                      % As the beamwidth is defined to be 4 standard deviations, i.e. 2 in each direction
            sigma_range_list    = sigma_range_0 + sigma_radial_list .* tan(incidence_angle_list);

            % Resulting range noise
            SN_noise_list       = normrnd(0, 1, [number_points, 1]);
            range_noise_list    = range_bias + sigma_range_list .* SN_noise_list;
            point_range_list_n  = point_range_list + range_noise_list;

            % Noisy point locations
            point_cloud_matrix_n = scanner_location + point_range_list_n .* vector_matrix;

            % Points exceeding the maximum incidence angle are removed
            incidence_angle_bool    = incidence_angle_list < max_incidence_angle;
            point_cloud_matrix_n    = point_cloud_matrix_n(incidence_angle_bool, :);
            point_cloud_cell_n{s}   = point_cloud_matrix_n;
        end

        % Noisy point cloud structure
        Number_Points_fun       = @(point_cloud_matrix) size(point_cloud_matrix, 1);
        number_points_list_n    = cellfun(Number_Points_fun, point_cloud_cell_n);

        Point_Cloud_Data_n = struct('point_cloud_cell', {point_cloud_cell_n}, 'number_points_list', number_points_list_n);

    %% Plot %%
        if Plot == true
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])  
    
            hold on
            grid on

            % The original point cloud
            point_cloud_matrix_total = vertcat(point_cloud_cell{:});
            scatter3(point_cloud_matrix_total(:, 1), point_cloud_matrix_total(:, 2), point_cloud_matrix_total(:, 3), 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', 'DisplayName', 'Original PC');

            point_cloud_matrix_total_n = vertcat(point_cloud_cell_n{:});
            scatter3(point_cloud_matrix_total_n(:, 1), point_cloud_matrix_total_n(:, 2), point_cloud_matrix_total_n(:, 3), 10, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'DisplayName', 'Noisy PC');

            % Axes
            xlabel('x [m]');
            ylabel('y [m]');
            zlabel('z [m]');
    
            axis equal
    
            view(45, 45);
    
            % Legend
            legend('show', 'location', 'eastoutside');

            % Formatting
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
    
            hold off

            % Pause message
            disp('The noisy point cloud has been determined. The figure closes and script continues upon a key-press.');
            pause();

            close(1);
        end
end