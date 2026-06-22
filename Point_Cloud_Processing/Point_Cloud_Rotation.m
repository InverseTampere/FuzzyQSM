% The point cloud, distributions and scanner locations are rotated by the given rotation matrix
% Structures that do not need to be rotated may be left as empty inputs

function [Point_Cloud_Coord_r, Point_Cloud_Distributions_r, point_cloud_matrix_r, Scanning_Parameters_r, Scanner_loc_cell_r, number_scanners] = Point_Cloud_Rotation(rotation_matrix, Point_Cloud_Coord, Point_Cloud_Distributions, Scanning_Parameters)

    %% Rotation function %%
        Rotation_fun = @(data_matrix) (rotation_matrix * data_matrix')';

    %% Point cloud centering %% 
        if ~isempty(Point_Cloud_Coord)
            % The data is rotated according to the rotation matrix
            point_cloud_cell        = Point_Cloud_Coord.point_cloud_cell;
            point_cloud_cell_r      = cellfun(Rotation_fun, point_cloud_cell, 'UniformOutput', false);

            point_cloud_matrix_r    = vertcat(point_cloud_cell_r{:});

            % And appended to a new structure
            Point_Cloud_Coord_r                     = Point_Cloud_Coord;
            Point_Cloud_Coord_r.point_cloud_cell    = point_cloud_cell_r;
        else
            % An empty coordinate output is returned otherwise
            Point_Cloud_Coord_r = [];
        end

    %% Distribution centering %%
        if ~isempty(Point_Cloud_Distributions)            
            % The expected values are rotated
            distribution_mu_cell        = Point_Cloud_Distributions.distribution_mu_cell;
            distribution_mu_cell_r      = cellfun(Rotation_fun, distribution_mu_cell, 'UniformOutput', false);

            point_cloud_matrix_r        = vertcat(distribution_mu_cell_r{:});

            % The axes are rotated  
            distribution_axes_cell      = Point_Cloud_Distributions.distribution_axes_cell;
            distribution_axes_cell_r    = cellfun(Rotation_fun, distribution_axes_cell, 'UniformOutput', false);

            % The covariance matrices are rotated
            distribution_Sigma_cell     = Point_Cloud_Distributions.distribution_Sigma_cell;
            Covariance_rotation_fun     = @(Sigma) rotation_matrix * Sigma * rotation_matrix';                                  % The covariance matrix is not simply expressed in x,y,z
            distribution_Sigma_cell_r   = cellfun(Covariance_rotation_fun, distribution_Sigma_cell, 'UniformOutput', false);

            % The new structure
            Point_Cloud_Distributions_r                         = Point_Cloud_Distributions;
            Point_Cloud_Distributions_r.distribution_mu_cell    = distribution_mu_cell_r;
            Point_Cloud_Distributions_r.distribution_axes_cell  = distribution_axes_cell_r;
            Point_Cloud_Distributions_r.distribution_Sigma_cell = distribution_Sigma_cell_r;

        else
            Point_Cloud_Distributions_r = [];
        end

    %% Centering the scanner locations %%
        if ~isempty(Scanning_Parameters)
            % Rotating the scanner locations
            Scanner_loc_cell        = Scanning_Parameters.Scanner_loc_cell;
            number_scanners         = Scanning_Parameters.number_scanners;
            
            Scanner_loc_cell_r      = cellfun(Rotation_fun, Scanner_loc_cell, 'UniformOutput', false);

            % The new structure
            Scanning_Parameters_r                   = Scanning_Parameters;
            Scanning_Parameters_r.Scanner_loc_cell  = Scanner_loc_cell_r;
        else
            [Scanning_Parameters_r, Scanner_loc_cell_r, number_scanners] = deal([]);
        end
end