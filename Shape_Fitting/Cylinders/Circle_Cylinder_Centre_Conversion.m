% This script converts between the 3D cylinder centre and the 2D circle centre, which is the projection on the cross-sectional plane with the point cloud centroid as origin

% Note: The circle centre and height are uniquely defined for each scanner location and given in each row of the output matrix
% Note: The circle centre (with corresponding height) may be given as an empty input

function [circle_centre_matrix, circle_height_list, cylinder_centre, cylinder_centre_c] = Circle_Cylinder_Centre_Conversion(circle_centre_matrix, circle_height_list, cylinder_centre, cylinder_direction, point_cloud_centroid, Scanner_loc_cell)

    %% Conversion vector basis %%
        % For consistency the point cloud centroid is used
        cylinder_vector_base_cell = Vector_Basis_Cylinder_Cross_Section_Projection(point_cloud_centroid, cylinder_direction, Scanner_loc_cell);

    %% Circle to cylinder centre %%
        if ~isempty(circle_centre_matrix)
            % The height above the cross-sectional plane is added in the third dimension
            circle_centre_matrix_3D = [circle_centre_matrix, circle_height_list];
            
            % The cylinder centre is computed using the first scanner's data
            cylinder_vector_basis   = cylinder_vector_base_cell{1};
            circle_centre_3D        = circle_centre_matrix_3D(1, :);

            cylinder_centre_c       = (cylinder_vector_basis' * circle_centre_3D')';
            cylinder_centre         = cylinder_centre_c + point_cloud_centroid;
        end

    %% Cylinder to circle centre %%
        if ~isempty(cylinder_centre)
            % Subtraction of the centroid
            cylinder_centre_c       = cylinder_centre - point_cloud_centroid;

            % Rotation by the vector bases
            num_dim                 = length(cylinder_centre);
            number_scanners         = length(Scanner_loc_cell);

            circle_centre_matrix    = zeros(number_scanners, num_dim - 1);
            circle_height_list      = zeros(number_scanners, 1);

            for s = 1 : number_scanners
                % This scanner's vector basis
                cylinder_vector_basis   = cylinder_vector_base_cell{s};
                circle_centre_3D        = (cylinder_vector_basis * cylinder_centre_c')';
    
                % Separation of the first two and third dimensions
                circle_centre_matrix(s, :)  = circle_centre_3D(1 : num_dim - 1);
                circle_height_list(s)       = circle_centre_3D(num_dim);
            end
        end

end