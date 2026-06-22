% This script creates coordinates in (x, y, z) of the given cylinder, discretised in circles spanning its length

function cylinder_coord_matrix = Cylinder_Coordinate_Generator(cyl_centre, cyl_radius, cyl_length, cyl_direction, number_circle_steps, number_length_steps)

    %% Cylinder coordinates %%
        % Unit vector length for the cylinder vector
        cyl_direction = cyl_direction / norm(cyl_direction);

        % The coordinates are generated for each circle along the length of the cylinder        
        length_steps        = cyl_length/2 * linspace(-1, 1, number_length_steps);
        cylinder_coord_cell = cell(1, number_length_steps);

        for l = 1 : number_length_steps
            % Centre of the circle for this length
            length          = length_steps(l);
            circle_centre   = cyl_centre + length*cyl_direction;

            % Circle coordinates at this length step
            circle_coord_matrix     = Circle_Coordinates(circle_centre, cyl_radius, cyl_direction, number_circle_steps);
            cylinder_coord_cell{l}  = circle_coord_matrix;
        end
    
        % The coordinates are merged
        cylinder_coord_matrix = vertcat(cylinder_coord_cell{:});
end