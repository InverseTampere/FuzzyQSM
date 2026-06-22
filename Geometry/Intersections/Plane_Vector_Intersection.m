% This script computes the intersection point between a plane and the given vectors
% For vectors parallel to the plane, the intersection and delta are defined as NaN

function [intersection_matrix, delta_list, incidence_angle_list, parallel_vectors] = Plane_Vector_Intersection(plane_point, plane_normal_vector, vector_start_matrix, vector_matrix)
    
    %% Manual inputs %%
        Plot    = false;     % [true, false] Shows the plane, vectors and their intersections

    %% Intersection points %%
        % The plane constant
        plane_constant = dot(plane_normal_vector, plane_point);
    
        % The distances along the vectors at which the intersections occur
        vector_dot_products = sum(plane_normal_vector .* vector_matrix, 2);
        start_dot_products  = sum(plane_normal_vector .* vector_start_matrix, 2);
        delta_list          = (plane_constant - start_dot_products) ./ vector_dot_products;
        
        % The incidence angle is the minimum angle between 0 and pi/2
        incidence_angle_list    = acos(abs(vector_dot_products));
        incidence_angle_list    = min(incidence_angle_list, pi - incidence_angle_list);

        % The intersection points
        intersection_matrix = vector_start_matrix + delta_list .* vector_matrix;

        % Note that when the dot product between the vectors is zero, the vectors are parallel to the plane and thus intersections will never take place
        % These are given a NaN value
        numerical_margin                            = 1e-3;
        parallel_vectors                            = abs(vector_dot_products) < numerical_margin;
        intersection_matrix(parallel_vectors, :)    = NaN;
        delta_list(parallel_vectors)                = NaN;

    %% Plot %%
        if Plot == true
            % Corner points of the plane            
            valid_intersection_matrix   = intersection_matrix(~parallel_vectors, :);
            plane_corner_matrix         = Plane_Corner_Points(plane_normal_vector, plane_point, valid_intersection_matrix);
        
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])     
            
            hold on
            grid on
            
            % The plane
            patch(plane_corner_matrix(:, 1), plane_corner_matrix(:, 2), plane_corner_matrix(:, 3), 'b', 'FaceAlpha', 0.5, 'DisplayName', 'Plane');
            
            % The normal vector to the plane
            plot3(plane_point(1) + [0, plane_normal_vector(1)], plane_point(2) + [0, plane_normal_vector(2)], plane_point(3) + [0, plane_normal_vector(3)], 'color', 'b', 'LineWidth', 2, 'DisplayName', 'Plane normal');
    
            % The intersection points
            scatter3(intersection_matrix(:, 1), intersection_matrix(:, 2), intersection_matrix(:, 3), 'filled', 'MarkerFaceColor' ,'r', 'DisplayName', 'Intersections');
            
            % The vectors
            number_vectors = length(delta_list);

            for v = 1 : number_vectors
                % The data of this vector
                vector          = vector_matrix(v, :); 
                vector_start    = vector_start_matrix(v, :);
                delta           = delta_list(v);
                
                pl = plot3(vector_start(1) + 2*delta*vector(1)*[-1, 1], vector_start(2) + 2*delta*vector(2)*[-1, 1], vector_start(3) + 2*delta*vector(3)*[-1, 1], 'color', 'r', 'LineWidth', 1, 'DisplayName', 'Vector');
                
                if v > 1
                    pl.HandleVisibility = 'Off';
                end
            end
            
            % Axes
            xlabel('x');
            ylabel('y');
            zlabel('z');
            
            axis equal
            
            view(45, 45);
            
            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    

            % Time to look at the figure
            disp('The figure is closed when a button is pressed, and the script continues.');
            pause();
            close(1);
        end
    
end