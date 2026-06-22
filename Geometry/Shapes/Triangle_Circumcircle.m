% The circumcircle of the given triangle is computed, for 2D or 3D
% triangle_vertices is expected to contain the coordinates of its three vertices in a 3x2 or 3x3 matrix (one vertex per row)
% If the vertices are (nearly) collinear, the returned values are NaN

function [circumcentre, circumradius] = Triangle_Circumcircle(triangle_vertices)

    %% Manual inputs %%
        numerical_margin    = 1e-3;         % [-] To check for collinearity
        Plot                = false;         % [true, false] To plot the triangle and its circumcircle

    %% Circumcircle %%
        % Triangle vertex coordinates
        [a, b, c]   = deal(triangle_vertices(1, :), triangle_vertices(2, :), triangle_vertices(3, :));
        num_dim     = length(a);

        % Check for collinearity
        ac_vector = a - c;
        bc_vector = b - c;
        ab_vector = a - b;

        if num_dim == 2
            % The area triangle per Heron's formula
            ac_length = norm(ac_vector);
            bc_length = norm(bc_vector);
            ab_length = norm(ab_vector);

            semi_perimeter  = (ac_length + bc_length + ab_length) / 2;
            area_triangle   = sqrt(semi_perimeter * (semi_perimeter - ac_length)*(semi_perimeter - bc_length)*(semi_perimeter - ab_length));

            % If it is below the margin
            if area_triangle < numerical_margin
                circumcentre = NaN(1, num_dim);
                circumradius = NaN;
                
                return
            end
            
        elseif num_dim == 3
            % Collinear edges are parallel, thus their cross product has a near zero norm
            cross_product = cross(ac_vector, bc_vector);        
            relative_norm = 2*norm(cross_product) / (norm(ac_vector) + norm(bc_vector));
    
            if relative_norm < numerical_margin
                circumcentre = NaN(1, num_dim);
                circumradius = NaN;
    
                return
            end
        else
            error('The triangle circumcircle can only be determined in 2D or 3D');
        end

        % Circumcircle
        if num_dim == 2
            % Two-dimensional circumcircle
            A   = det([a, 1; b, 1; c, 1]);
            B   = det([a, norm(a)^2; b, norm(b)^2; c, norm(c)^2]);
            s_x = 1/2*det([norm(a)^2, a(2), 1; norm(b)^2, b(2), 1; norm(c)^2, c(2), 1]);
            s_y = 1/2*det([a(1), norm(a)^2, 1; b(1), norm(b)^2, 1; c(1), norm(c)^2, 1]);
        
            circumcentre = [s_x, s_y] / A;
            circumradius = sqrt(B/A + (s_x^2 + s_y^2)/A^2);
    
            % Tiny triangles may have a very small complex radius due to rounding, which is converted to 0
            if imag(circumradius) ~= 0
                circumradius = 0;
            end

        elseif num_dim == 3
            % Three-dimensional circumcircle
            triangle_orth_vector    = cross(ac_vector, bc_vector);
            circumradius            = norm(ac_vector)*norm(bc_vector)*norm(ac_vector - bc_vector) / (2*norm(triangle_orth_vector));
            circumcentre            = cross(norm(ac_vector)^2*bc_vector - norm(bc_vector)^2*ac_vector, triangle_orth_vector) / (2*norm(triangle_orth_vector)^2) + c;
        end

    %% Plot %%
        if Plot == true
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])  
        
            hold on
            grid on
        
            % Triangle
            triangle_points = [a; b; c; a];     % To make it closed

            if num_dim == 2
                plot(triangle_points(:, 1), triangle_points(:, 2), 'LineWidth', 2, 'color', 'b');
            else
                plot3(triangle_points(:, 1), triangle_points(:, 2), triangle_points(:, 3), 'LineWidth', 2, 'color', 'b');
            end
        
            % Circumcircle
            number_coord = 1e2;

            if num_dim == 2
                circle_coord_matrix = Circle_Coordinates(circumcentre, circumradius, [], number_coord);
                circle_coord_matrix = [circle_coord_matrix; circle_coord_matrix(1, :)];                     % To make it closed
                plot(circle_coord_matrix(:, 1), circle_coord_matrix(:, 2), 'LineWidth', 2, 'color', 'r');
            else
                [unit_sphere_x_list, unit_sphere_y_list, unit_sphere_z_list] = sphere(number_coord);
                [sphere_x_list, sphere_y_list, sphere_z_list] = deal(circumradius*unit_sphere_x_list + circumcentre(1), circumradius*unit_sphere_y_list + circumcentre(2), circumradius*unit_sphere_z_list + circumcentre(3));

                surf(sphere_x_list, sphere_y_list, sphere_z_list, 'FaceColor', 'r', 'FaceAlpha', 0.25);
            end
        
            % Axes
            xlabel('x [m]');
            ylabel('y [m]');

            if num_dim == 3
                zlabel('z [m]');
                view(45, 45);
            end

            axis equal
        
            % Font
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
        
            hold off

            % Pause message
            disp('The circumcircle has been determined. The figure closes and script ends upon a key-press.');
            pause();

            close(1);
        end

end