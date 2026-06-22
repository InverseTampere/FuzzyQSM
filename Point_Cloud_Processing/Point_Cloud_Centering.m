% The point cloud and optionally distributions and scanner locations are given w.r.t. the point cloud centroid
% The scanner and alternatively point cloud coordinate / distribution inputs may be empty

function [Point_Cloud_Coord_c, Point_Cloud_Distributions_c, point_cloud_matrix_c, point_cloud_centroid, Scanning_Parameters_c, Scanner_loc_cell_c, number_scanners] = Point_Cloud_Centering(Point_Cloud_Coord, Point_Cloud_Distributions, Scanning_Parameters)

    %% Point cloud centering %% 
        if ~isempty(Point_Cloud_Coord)
            point_cloud_cell = Point_Cloud_Coord.point_cloud_cell;

            % Its centroid
            point_cloud_matrix      = vertcat(point_cloud_cell{:});
            point_cloud_centroid    = mean(point_cloud_matrix, 1);
    
            % Centering the point cloud
            point_cloud_matrix_c    = point_cloud_matrix - point_cloud_centroid;            
            Centering_fun           = @(matrix) matrix - point_cloud_centroid;
            point_cloud_cell_c      = cellfun(Centering_fun, point_cloud_cell, 'UniformOutput', false);
        
            % The new structure
            Point_Cloud_Coord_c                     = Point_Cloud_Coord;
            Point_Cloud_Coord_c.point_cloud_cell    = point_cloud_cell_c;

        else
            % An empty coordinate output is returned
            Point_Cloud_Coord_c     = [];
        end

    %% Distribution centering %%
        if ~isempty(Point_Cloud_Distributions)
            % The expected value of each distribution is used as the point cloud
            distribution_mu_cell    = Point_Cloud_Distributions.distribution_mu_cell;
            point_cloud_matrix      = vertcat(distribution_mu_cell{:});

            point_cloud_centroid    = mean(point_cloud_matrix, 1);
            point_cloud_matrix_c    = point_cloud_matrix - point_cloud_centroid;
            
            % The expected values are shifted by the centroid
            Centering_fun           = @(distribution_mu) distribution_mu - point_cloud_centroid;
            distribution_mu_cell_c  = cellfun(Centering_fun, distribution_mu_cell, 'UniformOutput', false);

            % The new structure
            Point_Cloud_Distributions_c                         = Point_Cloud_Distributions;
            Point_Cloud_Distributions_c.distribution_mu_cell    = distribution_mu_cell_c;
        else
            Point_Cloud_Distributions_c = [];
        end

    %% Centering the scanner locations %%
        if ~isempty(Scanning_Parameters)
            % The scanner locations
            Scanner_loc_cell    = Scanning_Parameters.Scanner_loc_cell;
            number_scanners     = Scanning_Parameters.number_scanners;
            
            % Centering
            Scanner_Centering_fun   = @(scanner_loc) scanner_loc - point_cloud_centroid;
            Scanner_loc_cell_c      = cellfun(Scanner_Centering_fun, Scanner_loc_cell, 'UniformOutput', false);

            % The new structure
            Scanning_Parameters_c                   = Scanning_Parameters;
            Scanning_Parameters_c.Scanner_loc_cell  = Scanner_loc_cell_c;
        else
            [Scanning_Parameters_c, Scanner_loc_cell_c, number_scanners] = deal([]);
        end
end