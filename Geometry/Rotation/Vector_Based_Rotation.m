% The given data is rotated around the rotation point, such that the vector forms the new z-axis
% The vector basis is returned, which includes two vectors that are parallel to the plane

function [rotated_points_matrix, vector_basis, plane_vectors] = Vector_Based_Rotation(points_matrix, vector, rotation_point)

    %% Vector basis %%
        % The vector is normalised
        vector  = vector / norm(vector);
        num_dim = length(vector);

        % If the normal vector is aligned with one of the axes, the vector basis is a (permuted, sign-shifted) identity matrix
        if max(abs(vector)) > 1 - 1e-6                                                        % Some margin due to computational rounding
            % The normal vector is put in the third dimension
            vector_basis                = zeros(num_dim);
            vector_basis(num_dim, :)    = vector;
                        
            % Aligned axis and first planar vector
            aligned_axis        = find(abs(vector) == max(abs(vector)));
            planar_axes         = setdiff(1 : num_dim, aligned_axis);
            planar_axis_1       = planar_axes(1);

            identity_matrix     = eye(num_dim);
            planar_vec_1        = identity_matrix(planar_axis_1, :);
            vector_basis(1, :)  = planar_vec_1;

            % The second planar vector is the cross product which can be positive or negative
            planar_vector_2     = Cross_Product_3D(planar_vec_1, planar_axis_1, vector, aligned_axis);
            vector_basis(2, :)  = planar_vector_2;

            plane_vectors = vector_basis(1:num_dim - 1, :);

        % Otherwise, two orthogonal vectors are found on the plane
        else
            plane_vector_a  = [-vector(2), vector(1), 0] / sqrt(vector(1)^2 + vector(2)^2);                  
            plane_vector_b  = cross(vector, plane_vector_a);       
            plane_vectors   = [plane_vector_a; plane_vector_b];

            vector_basis    = [plane_vectors; vector];
        end
    
    %% Point matrix rotation %%
        % The matrix is translated w.r.t. the rotation point
        points_matrix_transl    = points_matrix - rotation_point;
        
        % Then rotated using the vector basis
        rotated_points_matrix   = (vector_basis * points_matrix_transl')';
        
        % And translated back
        rotated_points_matrix   = rotated_points_matrix + rotation_point;

end