% This script generates the coordinates of a 3D ellipsoid
% Note that the returned number of coordinates may slightly differ from the given number

% The coordinate frame [a, b, c] is aligned with [x, y, z] before rotation
% The radii are expected to be given in order [a, b, c] NOT in order of magnitude
% The ellipsoid axes are row vectors of the rotated a, b, c vectors
% Note that the rotation matrix is then simply the inverse of the ellipsoid axes

function [ellipsoid_coord_matrix, number_coords] = Ellipsoid_Coordinate_Generator(ellipsoid_centre, ellipsoid_radii, ellipsoid_axes, number_coords)

    %% Unit sphere %%
        % To create the ellipsoid, a unit sphere at the origin is created first
        n_ellipse = round(sqrt(number_coords));                                 % As it returns an n*n matrix, the square root is taken
        [x_unit_sphere, y_unit_sphere, z_unit_sphere] = sphere(n_ellipse);
        
        unit_sphere_coord_matrix    = [x_unit_sphere(:), y_unit_sphere(:), z_unit_sphere(:)];
        number_coords               = size(unit_sphere_coord_matrix, 1);        % Because of the earlier rounding, this may differ from the given number
        
    %% Ellipsoid %%
        % The ellipsoid is formed by multiplying the coordinates with the given radii and moving it to its centre
        ellipsoid_coord_matrix_e    = ellipsoid_radii .* unit_sphere_coord_matrix;

        % The final ellipsoid is then attained through rotation and trnaslation
        ellipsoid_coord_matrix_t    = (ellipsoid_axes \ ellipsoid_coord_matrix_e')';
        ellipsoid_coord_matrix      = ellipsoid_coord_matrix_t + ellipsoid_centre;

end