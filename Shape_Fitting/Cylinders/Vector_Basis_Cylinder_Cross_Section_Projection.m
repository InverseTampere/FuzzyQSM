% Vector bases are computed for the given cylinder and scanner locations
% The vector bases are such that the z-axis aligns with the cylinder axis and the scanner lies on the y-z plane, i.e.
% the x and y axes form the cross-sectional plane with the propagation axis in the -y direction

function cylinder_vector_basis_cell = Vector_Basis_Cylinder_Cross_Section_Projection(cylinder_centre, cylinder_direction, Scanner_loc_cell)

    %% Vector bases %%
        % Normalising the input vector
        cylinder_direction = cylinder_direction / norm(cylinder_direction);

        % Vector bases for each scanner
        number_scanners             = length(Scanner_loc_cell);
        cylinder_vector_basis_cell  = cell(1, number_scanners);

        for s = 1 : number_scanners
            % Relative scanner location w.r.t. the cylinder centre
            scanner_loc_c = Scanner_loc_cell{s} - cylinder_centre;

            % Scanner vector basis
            a_vector = cross(cylinder_direction, scanner_loc_c);
            a_vector = a_vector / norm(a_vector);
    
            b_vector = cross(cylinder_direction, a_vector);
    
            if dot(b_vector, scanner_loc_c) < 0             % To ensure that the y-coordinate of the scanner is positive
                a_vector = -a_vector;
                b_vector = -b_vector;
            end
    
            vector_basis                    = [a_vector; b_vector; cylinder_direction];
            cylinder_vector_basis_cell{s}   = vector_basis;
        end
end