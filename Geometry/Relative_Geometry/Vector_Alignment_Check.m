% This script checks whether the alignment between a given direction vector and a set of vectors exceeds the given threshold in radians
% Note that it returns normalised vectors

function [aligned_vectors, oblique_vectors, aligned_boolean, theta_list] = Vector_Alignment_Check(direction_vector, vector_matrix, alignment_threshold)

    %% Vectors are normalised %%
        % This ensures that no division is necessary when computing the angle
        direction_vector    = direction_vector / norm(direction_vector);
        vector_matrix       = vector_matrix ./ sqrt(sum(vector_matrix.^2, 2));

    %% Angle between the given vector and set of vectors %%
        % The direction vector is repeated for the dot product
        num_vectors             = size(vector_matrix, 1);
        direction_vector_matrix = repmat(direction_vector, [num_vectors, 1]);

        % The angle is irrespective of the orientation of the given vector
        theta_list_plus     = acos(dot(vector_matrix, direction_vector_matrix, 2));     
        theta_list_min      = acos(dot(vector_matrix, -direction_vector_matrix, 2));

        theta_list          = min(abs(theta_list_plus), abs(theta_list_min));

    %% Set of aligned and oblique vectors %%
        % Aligned vectors are those whose angle to the given vector is too small
        aligned_boolean     = (theta_list < alignment_threshold)';
        
        aligned_vectors     = vector_matrix(aligned_boolean, :);
        oblique_vectors     = vector_matrix(~aligned_boolean, :);
    
end