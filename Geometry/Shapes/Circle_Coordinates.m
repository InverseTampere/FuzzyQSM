% This script computes the coordinates of a circle, given the radius and centre location.
% If the normal vector is given, 3D coordinates are returend. Otherwise 2D.
% Note that they start at the top of the circle, and run clockwise.

function circle_coord_matrix = Circle_Coordinates(circle_centre, circle_radius, circle_normal_vector, number_coord)

    %% Circle coordinates %%
        % Angles at the desired discretisation level
        theta_list                      = linspace(0, 2*pi, number_coord + 1)';
        theta_list(number_coord + 1)    = [];       % To avoid duplicate points at 0 and 2pi, the final point is removed
        
        % The x and y coordinates, centered on the circle
        circle_coord_matrix_c = circle_radius * [cos(theta_list), sin(theta_list)];

        % Rotation, if applicable
        if ~isempty(circle_normal_vector)
            % A third dimension is added
            circle_coord_matrix_c   = [circle_coord_matrix_c, zeros(number_coord, 1)];

            % Rotation
            [~, vector_basis, ~]    = Vector_Based_Rotation(circle_centre, circle_normal_vector, circle_centre); 
            circle_coord_matrix_c   = (vector_basis \ circle_coord_matrix_c')';
        end

        % The centre is added
        circle_coord_matrix = circle_coord_matrix_c + circle_centre;

end