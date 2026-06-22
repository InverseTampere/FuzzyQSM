% This script retrieves the mesh and circles in the .obj file produced by Markku Åkerblom's Blender add-on
% Note that it expects a singular object, i.e. for branch separation to NOT have been ticked in Blender

function [Smoothed_Model_Circle_Geometry, Smoothed_Model_Mesh] = Smoothed_Model_Circles(object_file_path)

    %% Inputs %%
        % Object file specifications
        delimiter               = ' ';                          % Delimiter between the object file entries, i.e. ' '
        Y_Axis_Up               = true;                        % [true, false] Blender permutes the axes by default
        Plot                    = false;                        % [true, false] Shows the object geometry

        % Circle geometry
        number_circle_points    = 20;                           % [-] Number of points per circle
        allowed_circle_error    = 1.0;                          % [%] Allowed relative error in the circle radius
        convergence_threshold   = 1e-3;                         % [-] Convergence threshold for finding the circle normal vectors

    %% Read the data %%
        % The mesh is recovered from the object file
        [Vertices, Triangles]   = Object_File_Reader(object_file_path, delimiter, Y_Axis_Up, Plot);
        Smoothed_Model_Mesh     = struct('Vertices', Vertices, 'Triangles', Triangles);

        vertex_matrix               = Vertices.coordinates;
        [number_vertices, num_dim]  = size(vertex_matrix);

    %% Determine circle geometry %%
        number_circles              = ceil(number_vertices / number_circle_points);
    
        circle_radius_list          = zeros(number_circles, 1);
        circle_centre_matrix        = zeros(number_circles, num_dim);
        circle_normal_vector_matrix = zeros(number_circles, num_dim);
        circle_start_ind_list       = zeros(number_circles, 1);
        circle_end_ind_list         = zeros(number_circles, 1);

        end_ind = 0;

        for c = 1 : number_circles
            % This circle's points
            start_ind       = end_ind + 1;
            end_ind         = start_ind + number_circle_points - 1;
            circle_points   = vertex_matrix(start_ind : end_ind, :);
    
            % Centre
            circle_centre               = mean(circle_points, 1);
            circle_centre_matrix(c, :)  = circle_centre;
    
            % Radius
            norm_list               = sqrt(sum((circle_points - circle_centre).^2, 2));
            circle_radius           = mean(norm_list);
            circle_radius_list(c)   = circle_radius;
    
            circle_error = max(abs(norm_list - circle_radius)) / circle_radius * 100;               % Maximum relative error from this radius
            
            numerical_margin = 1e-3;                                                                % Tiny circles can deviate more
            if circle_radius > numerical_margin && circle_error > allowed_circle_error
                % It is possible for the circle to be incomplete
                diff_number_points = number_circles*number_circle_points - number_vertices;
                if diff_number_points > 0
                    % Adjusted end point
                    end_ind         = end_ind - diff_number_points;
                    circle_points   = vertex_matrix(start_ind : end_ind, :);
    
                    % The centre and radius now have to be found using least-squares
                    Sigma                   = cov(circle_points);
                    [~, normal_vector]      = Power_Iteration(Sigma, convergence_threshold);

                    origin                              = zeros(1, num_dim); 
                    [circle_points_r, vector_basis, ~]  = Vector_Based_Rotation(circle_points, normal_vector, origin);
                    circle_points_2D                    = circle_points_r(:, 1 : num_dim - 1);
                    height                              = circle_points_r(1, num_dim);
                    
                    weight_list                         = ones(end_ind - start_ind + 1, 1);
                    Method                              = 'Algebraic';
                    LS_Diagnostics                      = false;
                    [circle_centre_2D, circle_radius]   = Least_Squares_Circle_Fitting(circle_points_2D, weight_list, Method, LS_Diagnostics);
                    
                    circle_centre_r = [circle_centre_2D, height];
                    circle_centre   = (vector_basis' * circle_centre_r')';

                    norm_list = sqrt(sum((circle_points - circle_centre).^2, 2));

                    circle_error = max(abs(norm_list - circle_radius)) / circle_radius * 100;               % Maximum relative error from this radius

                    if circle_radius > numerical_margin && circle_error > allowed_circle_error
                        error('The circle error is %.3g %%.', circle_error);
                    end
                else
                    error('The circle error is %.3g %%.', circle_error);
                end
            end
    
            % Saving the start and end indices as they correspond to the 
            circle_start_ind_list(c)    = start_ind;
            circle_end_ind_list(c)      = end_ind;

            % Normal vector
            Sigma               = cov(circle_points);
            [~, normal_vector]  = Power_Iteration(Sigma, convergence_threshold);

            circle_normal_vector_matrix(c, :) = normal_vector;
        end
       
        % Create the structure
        Smoothed_Model_Circle_Geometry = struct('number_circle_points', number_circle_points, 'radius', circle_radius_list, 'centre', circle_centre_matrix, 'normal_vector', circle_normal_vector_matrix, 'number_circles', number_circles, 'vertex_start_indices', circle_start_ind_list, 'vertex_end_indices', circle_end_ind_list);
end
