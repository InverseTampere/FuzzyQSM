% This script creates samples that are equally spaced and covering a circular area in 2D or 3D (requiring a normal vector)
% The total number of samples may differ slightly from the initially given value

function [sample_coord_matrix, number_samples, delta_area] = Equal_Area_Circular_Sampler(circle_radius, circle_centre, normal_vector, number_samples)
        
    %% Circle on x-y plane %%
        % The number of samples in each circle is equivalent to a square, with the following samples per side of each quarter
        n_quarter = floor(2/pi * sqrt(number_samples));

        % The central row and column
        var_list    = circle_radius/(n_quarter + 1) * linspace(1, n_quarter + 1, n_quarter);
        
        x_list_c    = [-var_list, 0, var_list];         % Note how the origin is added as well
        y_list_c    = [-var_list, 0, var_list];
        
        % The full grid
        [x_matrix, y_matrix] = meshgrid(x_list_c, y_list_c);
        
        % Of all the samples, only the ones with norm within the radius are kept
        norm_matrix     = sqrt(x_matrix.^2 + y_matrix.^2);
        
        x_circle_list   = x_matrix(norm_matrix <= circle_radius);
        y_circle_list   = y_matrix(norm_matrix <= circle_radius);
        
        number_samples  = length(x_circle_list);
        
        % The area per sample is then
        delta_area = pi * circle_radius^2 / number_samples;
        
    %% Rotation and translation %%
        % Translated according to the centroid
        sample_coord_matrix = [x_circle_list, y_circle_list];

        % If there is no given normal vector the data is presumed 2D
        if isempty(normal_vector)
            sample_coord_matrix = sample_coord_matrix + circle_centre;
            
        % Otherwise the samples are rotated s.t. the normal vector forms the z-axis        
        else
            % A third dimension is appended
            z_circle_list           = zeros(number_samples, 1);
            sample_coord_matrix_r   = [sample_coord_matrix, z_circle_list];
            
            % They are rotated
            [~, vector_basis, ~]    = Vector_Based_Rotation(sample_coord_matrix_r, normal_vector, circle_centre);
            sample_coord_matrix_t   = (vector_basis' * sample_coord_matrix_r')';
            sample_coord_matrix     = sample_coord_matrix_t + circle_centre;
        end
end