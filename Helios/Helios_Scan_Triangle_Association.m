% If the Helios scan was performed without noise, this script determines which triangle each point intersected

function triangle_index_cell = Helios_Scan_Triangle_Association(point_cloud_cell, Scanning_Parameters, Triangular_Mesh)
    
    %% Structure inputs %%
        % Scanning parameters
        Scanner_loc_cell                = Scanning_Parameters.Scanner_loc_cell;
        number_scanners                 = Scanning_Parameters.number_scanners;

        % Triangular mesh
        vertex_coord_matrix             = Triangular_Mesh.Vertices.coordinates;
        triangle_vertex_index_matrix    = Triangular_Mesh.Triangles.vertex_indices;

        clearvars Triangular_Mesh

    %% Manual inputs %%
        Diagnostics                     = false;        % [true, false] Can be used to check the code that associates an individual point with the triangles
    
    %% Determining each point's triangle %%
        % The range margin of each triangle follows from half its maximum vertex length
        range_margin_list = Triangle_Range_Margin(vertex_coord_matrix, triangle_vertex_index_matrix);

        % Each scanner is dealt with separately
        triangle_index_cell = cell(1, number_scanners);

        for s = 1 : number_scanners
            % This scanner's data
            scanner_location    = Scanner_loc_cell{s};
            point_cloud_matrix  = point_cloud_cell{s};

            % Vectors to the points and their ranges
            vector_matrix       = point_cloud_matrix - scanner_location;
            point_range_list    = sqrt(sum(vector_matrix.^2, 2));

            % Ranges to the vertices
            vertex_vector_matrix    = vertex_coord_matrix - scanner_location;
            vertex_range_list       = sqrt(sum(vertex_vector_matrix.^2, 2));
            
            % Triangle ranges - including the radius
            triangle_range_matrix   = vertex_range_list(triangle_vertex_index_matrix);
            triangle_range_LB_list  = min(triangle_range_matrix, [], 2) - range_margin_list;
            triangle_range_UB_list  = max(triangle_range_matrix, [], 2) + range_margin_list;

            clearvars vector_matrix vertex_vector_matrix vertex_range_list triangle_range_matrix
        
            % Points can only lie within triangles that cover their range
            Triangle_Candidates_fun     = @(point_range) find(point_range >= triangle_range_LB_list & point_range <= triangle_range_UB_list);
        
            number_points               = length(point_range_list);
            triangle_candidates_cell    = cell(number_points, 1);

            for p = 1 : number_points
                point_range                 = point_range_list(p);
                triangle_candidates         = Triangle_Candidates_fun(point_range);
                triangle_candidates_cell{p} = triangle_candidates;
            end

            % If no candidates were found for a triangle, there is a problem
            number_candidates_list = cellfun(@length, triangle_candidates_cell);
        
            if min(number_candidates_list) == 0
                % Vectors without any triangle candidates
                zero_candidates_bool = number_candidates_list == 0;
    
                %--% Figure for diagnostic purposes %--%
                figure(1)
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])  
        
                hold on
                grid on
        
                sgtitle('Vectors w/o triangle candidates');
        
                % The triangular mesh
                trimesh(triangle_vertex_index_matrix, vertex_coord_matrix(:, 1), vertex_coord_matrix(:, 2), vertex_coord_matrix(:, 3), 'LineWidth', 1, 'FaceColor', 'none', 'EdgeColor', 'k');
        
                % Vectors w/o candidates
                for p = 1 : number_points
                    if zero_candidates_bool(p) == true
                        point = point_cloud_matrix(p, :);
                        plot3([scanner_location(1), point(1)], [scanner_location(2), point(2)], [scanner_location(3), point(3)], 'LineWidth', 1, 'color', 'r');
                    end
                end
        
                % Points w/o candidates
                scatter3(point_cloud_matrix(zero_candidates_bool, 1), point_cloud_matrix(zero_candidates_bool, 2), point_cloud_matrix(zero_candidates_bool, 3), 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
        
                % Axes
                xlabel('x [m]');
                ylabel('y [m]');
                zlabel('z [m]');
        
                axis equal
        
                view(45, 45);
        
                % Formatting
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
        
                hold off
        
                %--% Error message %--%
                error('No triangles were found for %i intersects.', sum(zero_candidates_bool));
            end
            
            % Determine the correct triangle
            [number_points, num_dim]        = size(point_cloud_matrix);
            point_cell                      = mat2cell(point_cloud_matrix, ones(1, number_points), num_dim);

            Triangle_Finder_fun                 = @(point, triangle_candidate_bool) Triangle_Point_Association(point, vertex_coord_matrix, triangle_vertex_index_matrix(triangle_candidate_bool, :), Diagnostics);
            [found_triangle_candidate_cell, ~]  = cellfun(Triangle_Finder_fun, point_cell, triangle_candidates_cell, 'UniformOutput', false);
        
            Triangle_Index_fun      = @(triangle_candidates, candidate_index) triangle_candidates(candidate_index);
            triangle_index_list     = cellfun(Triangle_Index_fun, triangle_candidates_cell, found_triangle_candidate_cell);
            triangle_index_cell{s}  = triangle_index_list;
        end

    %% Local functions %%
        % Each triangle's range margin
        function range_margin_list = Triangle_Range_Margin(vertex_coord_matrix, triangle_vertex_index_matrix)
            % The triangle vertex indices
            [vertex_1_ind_list, vertex_2_ind_list, vertex_3_ind_list] = Column_Deal(triangle_vertex_index_matrix);

            % The edge lengths
            edge_12_vector_matrix = vertex_coord_matrix(vertex_1_ind_list, :) - vertex_coord_matrix(vertex_2_ind_list, :);
            edge_23_vector_matrix = vertex_coord_matrix(vertex_2_ind_list, :) - vertex_coord_matrix(vertex_3_ind_list, :);
            edge_13_vector_matrix = vertex_coord_matrix(vertex_1_ind_list, :) - vertex_coord_matrix(vertex_3_ind_list, :);

            edge_12_length_list = sqrt(sum(edge_12_vector_matrix.^2, 2));
            edge_23_length_list = sqrt(sum(edge_23_vector_matrix.^2, 2));
            edge_13_length_list = sqrt(sum(edge_13_vector_matrix.^2, 2));

            % Half the maximum length is used as the range margin
            max_edge_length_list    = max([edge_12_length_list, edge_23_length_list, edge_13_length_list], [], 2);
            range_margin_list       = max_edge_length_list / 2;
        end

        % Determine a point's triangle 
        function [triangle_index, min_uvw] = Triangle_Point_Association(point, vertex_coord_matrix, triangle_vertex_index_matrix, Diagnostics)
            % Points may be off the surface due to rounding errors
            distance_margin     = 1e-3;         % In metres
            barycentric_margin  = 0.50;         % Relative to the triangle size. Extremely flat triangles can have relatively large errors
    
            % Triangle vertices
            [a_ind_list, b_ind_list, c_ind_list]    = Column_Deal(triangle_vertex_index_matrix);
            [a_matrix, b_matrix, c_matrix]          = deal(vertex_coord_matrix(a_ind_list, :), vertex_coord_matrix(b_ind_list, :), vertex_coord_matrix(c_ind_list, :));
        
            % Normal vectors
            normal_vector_matrix    = cross(b_matrix - a_matrix, c_matrix - a_matrix, 2);
            parallelogram_area_list = sqrt(sum(normal_vector_matrix.^2, 2));
            normal_vector_matrix    = normal_vector_matrix ./ parallelogram_area_list;
    
            % Use the distance from the point normal to the plane to see which planes the point is near
            normal_distance_list    = dot(normal_vector_matrix, point - a_matrix, 2);           % Note that it is a signed distance
            near_plane_bool         = abs(normal_distance_list) < distance_margin;
    
            if sum(near_plane_bool) == 0
                % If there is no near plane, there is a problem
                End_Message_fun = @() error('No plane near to the projected point was found. The minimum distance is %.3g m', min(abs(normal_distance_list)));
            
                Diagnostic_Plotter(point, vertex_coord_matrix, triangle_vertex_index_matrix, [], End_Message_fun)

            elseif sum(near_plane_bool) == 1
                % If one candidate is found, that is selected as the point's triangle
                triangle_index  = find(near_plane_bool);
                min_uvw         = 1;
    
            else
                % Restrict the search to triangles whose plane the point is near
                a_matrix                = a_matrix(near_plane_bool, :);
                b_matrix                = b_matrix(near_plane_bool, :);
                c_matrix                = c_matrix(near_plane_bool, :);
                parallelogram_area_list = parallelogram_area_list(near_plane_bool);
                normal_vector_matrix    = normal_vector_matrix(near_plane_bool, :);
                normal_distance_list    = normal_distance_list(near_plane_bool);
    
                % Point projected onto these near planes
                proj_point_matrix       = point - normal_distance_list .* normal_vector_matrix;
    
                % Barycentric coordinates                
                u_list      = sqrt(sum(cross(c_matrix - b_matrix, proj_point_matrix - b_matrix, 2).^2, 2)) ./ parallelogram_area_list;
                v_list      = sqrt(sum(cross(a_matrix - c_matrix, proj_point_matrix - c_matrix, 2).^2, 2)) ./ parallelogram_area_list;
                w_list      = sqrt(sum(cross(b_matrix - a_matrix, proj_point_matrix - a_matrix, 2).^2, 2)) ./ parallelogram_area_list;

                % The correct triangle should have a sum of 1, and thus be the minimum
                [min_uvw, inside_index] = min(u_list + v_list + w_list);
    
                if min_uvw > 1 + barycentric_margin
                    % It may exceed it somewhat due to rounding, especially for very flat triangles, but it should not be excessive
                    End_Message_fun = @() error('The point was not found to be inside any triangle. The min. u+v+w is %.3g', min_uvw);

                    Diagnostic_Plotter(point, vertex_coord_matrix, triangle_vertex_index_matrix(near_plane_bool, :), [], End_Message_fun);    
                end

                % Converted to the original indices
                near_plane_triangles    = find(near_plane_bool);
                triangle_index          = near_plane_triangles(inside_index);
            end     

            % Diagnostic plot
            if Diagnostics == true
                Message_fun     = @() disp('The triangle has been found. The figure will close and script continue upon a key-press.');
                Pause_fun       = @() pause();
                Close_fun       = @() close(1);
                Fun_Cell        = {Message_fun, Pause_fun, Close_fun};
                End_Message_fun = @() cellfun(@feval, Fun_Cell);

                Diagnostic_Plotter(point, vertex_coord_matrix, triangle_vertex_index_matrix, triangle_index, End_Message_fun);
            end
        end

        % Function that shows the point and triangles
        function Diagnostic_Plotter(point, vertex_coord_matrix, triangle_vertex_index_matrix, triangle_index, End_Message_fun)

            %--% Plot %--%
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])  
    
            hold on
            grid on
            
            % Full mesh
            trimesh(triangle_vertex_index_matrix, vertex_coord_matrix(:, 1), vertex_coord_matrix(:, 2), vertex_coord_matrix(:, 3), 'LineWidth', 1, 'FaceColor', 'none', 'EdgeColor', 'k', 'DisplayName', 'Full mesh');
        
            % The corresponding triangle (if found)
            if ~isempty(triangle_index)
                trimesh(triangle_vertex_index_matrix(triangle_index, :), vertex_coord_matrix(:, 1), vertex_coord_matrix(:, 2), vertex_coord_matrix(:, 3), 'LineWidth', 1, 'FaceColor', 'none', 'EdgeColor', 'c', 'DisplayName', 'Corresponding triangle');
            end

            % Point
            scatter3(point(1), point(2), point(3), 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'DisplayName', 'Point');
    
            % Axes
            xlabel('x [m]');
            ylabel('y [m]');
            zlabel('z [m]');
    
            axis equal
    
            view(45, 45);
    
            % Legend
            legend('show', 'location', 'eastoutside');

            % Formatting
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
    
            hold off

            %--% End message %--%
            % Can be a simple pause or an error
            End_Message_fun();
        end
end