% This script provides a coordinate matrix that is rotated such that the given vector forms the new z-axis
% Note that the coordinates are expected to be given in column format, i.e. [X, Y, Z]

function [coordinate_matrix_rot, rotation_matrix] = Rotation_3D(coordinate_matrix, vector_dir, vector_start)
    
    %% Rotation matrix %%
        % If the given vector is parallel to the z-axis, rotation is not needed
        vector_dir_n    = vector_dir / norm(vector_dir);        % Normalised vector direction
        z_component     = vector_dir_n(3);                      % The component in the z-axis
            
        if abs(z_component) > 1 - 1e-6          % If the magnitude of the z-component is near 1, it is parallel to the z-axis
            Parallel_Vectors = true;
        else
            Parallel_Vectors = false;
        end
        
        % If the vectors are parallel, the rotation matrix is simply the identity matrix
        [~, num_dim] = size(coordinate_matrix);
                     
        if Parallel_Vectors == true
            rotation_matrix = eye(num_dim);
            
        else
            % The dot- and cross-products between the z-axis and vector direction
            normal_vector   = cross([0, 0, 1], vector_dir_n);
            vector_angle    = dot([0, 0, 1], vector_dir_n);

            % The skew-symmetric cross-product matrix of the normal vector
            U = [0,                 -normal_vector(3),  normal_vector(2);
                 normal_vector(3),  0,                  -normal_vector(1);
                 -normal_vector(2), normal_vector(1),   0];

             % The rotation matrix
             rotation_matrix = eye(num_dim) + U + U^2 / (1 + vector_angle);
             rotation_matrix = inv(rotation_matrix);
        end
        
     %% Coordinate matrix rotation %%
         % The coordinates are translated according to the start point of the vector
        coordinate_matrix_transl = coordinate_matrix - vector_start;

        % They are then rotated
        coordinate_matrix_rot   = (rotation_matrix * coordinate_matrix_transl')';
        
        % And translated according to the vector start
        coordinate_matrix_rot   = coordinate_matrix_rot + vector_start;

end