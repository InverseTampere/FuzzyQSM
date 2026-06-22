% This script generates the coordinates for surf to create a cylinder surface with

function [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, number_coord] = Cylinder_Surface_Generator(cyl_radius, cyl_length, cyl_centre, cyl_dir, number_coord)

    %% Cylinder coords %%
        % The number of steps is changed, as the coordinates are doubled later
        number_steps = floor(number_coord / 2);

        % Straight cylinder coordinates
        [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z] = cylinder(cyl_radius, number_steps);
        cylinder_coord_z = cyl_length * cylinder_coord_z - cyl_length/2;

        cylinder_coord_matrix = [cylinder_coord_x(:), cylinder_coord_y(:), cylinder_coord_z(:)];
        
        % Rotated
        [~, vector_basis]           = Rotation_3D(cylinder_coord_matrix, cyl_dir, cyl_centre);
        cylinder_coord_matrix_rot   = (vector_basis \ cylinder_coord_matrix')';
        
        % Translated
        cylinder_coord_matrix_trans = cylinder_coord_matrix_rot + cyl_centre;
        
        % Right matrix shape
        cylinder_coord_x = reshape(cylinder_coord_matrix_trans(:, 1), [2, number_steps + 1]);
        cylinder_coord_y = reshape(cylinder_coord_matrix_trans(:, 2), [2, number_steps + 1]);
        cylinder_coord_z = reshape(cylinder_coord_matrix_trans(:, 3), [2, number_steps + 1]);
        
        number_coord = length(cylinder_coord_x(:));
end
            