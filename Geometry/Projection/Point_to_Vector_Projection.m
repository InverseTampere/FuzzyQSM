% This script projects points onto a vector, regardless of dimension
% Delta is the distance from vector start to the projections along the vector
% Omega is the distance from the original points to their respective projections

function [projected_point_matrix, delta_list, omega_list] = Point_to_Vector_Projection(point_matrix, vector, vector_start)

    %% Point to vector projection %%
        % The vector from the start point of the vector to the projected points
        num_points      = size(point_matrix, 1);
        vector_matrix   = repmat(vector, [num_points, 1]);
        
        delta_list                  = dot(point_matrix - vector_start, vector_matrix, 2) / dot(vector, vector);        % The distance along the vector
        projection_vector_matrix    = delta_list .* vector;

        % The projected point is then translated by said vector's start point
        projected_point_matrix      = projection_vector_matrix + vector_start;
        
        % Distance from the original points to their projections
        omega_list = sqrt(sum((projected_point_matrix - point_matrix).^2, 2));
        
end