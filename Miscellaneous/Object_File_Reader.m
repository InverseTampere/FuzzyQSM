% This file reads Wavefront .obj files that are expected to contain only a single object
% Blender by default permutes the axes of the resulting data, causing the y-axis to be facing up and the third coordinate to be negative. 
% If this is the case, give Y_Axis_Up as true

function [Vertices, Faces] = Object_File_Reader(obj_file_path, delimiter, Y_Axis_Up, Plot)

    %% Vertex data %%
        % Read the file's lines
        object_file_lines = readlines(obj_file_path, 'EmptyLineRule', 'Skip');
    
        % Line splitting function
        Line_Splitter_fun = @(line) cellstr(strsplit(line, delimiter));                         % Note that it also converts the string array to cell array
    
        %--% Vertex coordinates %--%
        % Vertex lines
        vertex_line_bool    = startsWith(object_file_lines, sprintf('v%s', delimiter));         % The delimiter means vt and vn are excluded
        vertex_lines        = object_file_lines(vertex_line_bool);
        
        % The parts of each line
        vertex_line_parts   = arrayfun(Line_Splitter_fun, vertex_lines, 'UniformOutput', false);
        vertex_line_parts   = vertcat(vertex_line_parts{:});
    
        % The coordinates
        vertex_coord_cell   = vertex_line_parts(:, 2:4);
        vertex_coord_matrix = cellfun(@str2double, vertex_coord_cell);
    
        if Y_Axis_Up == true
            % Permutation to ensure the z-axis is the third coordinate and is pointing in the right direction
            vertex_coord_matrix         = vertex_coord_matrix(:, [1, 3, 2]);
            vertex_coord_matrix(:, 2)   = -vertex_coord_matrix(:, 2);
        end
    
        % If there were 5 parts to each line, they need to be scaled
        if size(vertex_line_parts, 2) == 5
            vertex_scale_cell   = vertex_line_parts(:, 5);
            vertex_scale_list   = cellfun(@str2double, vertex_scale_cell);
            vertex_coord_matrix = vertex_coord_matrix ./ vertex_scale_list;
        end
    
        %--% Texture coordinates %--%
        % Texture lines
        texture_line_bool   = startsWith(object_file_lines, 'vt');      
        texture_lines       = object_file_lines(texture_line_bool);
    
        if ~isempty(texture_lines)              % They are optional
            Textures            = true;
    
            % The parts of each line
            texture_line_parts  = arrayfun(Line_Splitter_fun, texture_lines, 'UniformOutput', false);
            texture_line_parts  = vertcat(texture_line_parts{:});
        
            % Texture coordinates may be 2D or 3D, so all columns after the first are selected
            texture_coord_cell      = texture_line_parts(:, 2 : end);
            texture_coord_matrix    = cellfun(@str2double, texture_coord_cell);
        else
            Textures                = false;
            texture_coord_matrix    = [];
        end
    
        %--% Vertex normals %--%
        % Vertex normal lines
        vertex_normal_line_bool = startsWith(object_file_lines, 'vn');      
        vertex_normal_lines     = object_file_lines(vertex_normal_line_bool);
    
        % The parts of each line
        vertex_normal_line_parts = arrayfun(Line_Splitter_fun, vertex_normal_lines, 'UniformOutput', false);
        vertex_normal_line_parts = vertcat(vertex_normal_line_parts{:});
    
        % The vector components
        vertex_normal_cell   = vertex_normal_line_parts(:, 2:4);
        vertex_normal_matrix = cellfun(@str2double, vertex_normal_cell);
    
        if Y_Axis_Up == true
            % Permutation to ensure the z-axis is the third coordinate
            vertex_normal_matrix = vertex_normal_matrix(:, [1, 3, 2]);
        end
    
        % Ensure they are unit vectors
        vector_norm_list        = sqrt(sum(vertex_normal_matrix.^2, 2));
        vertex_normal_matrix    = vertex_normal_matrix ./ vector_norm_list;
    
        %--% Vertex structure %--%
        Vertices = struct('coordinates', vertex_coord_matrix, 'textures', texture_coord_matrix, 'normal_vectors', vertex_normal_matrix);

    %% Face elements %%
        % Face element lines
        face_line_bool  = startsWith(object_file_lines, 'f');      
        face_lines      = object_file_lines(face_line_bool);
    
        % Vertex parts in each line
        face_vertex_parts = arrayfun(Line_Splitter_fun, face_lines, 'UniformOutput', false);
        face_vertex_parts = vertcat(face_vertex_parts{:});
        face_vertex_parts = face_vertex_parts(:, 2:4);              % The identifying f column is removed
    
        % The indices are split by /
        Face_Index_Splitter_fun = @(vertex_part) str2double(strsplit(vertex_part, '/'));
        face_vertex_index_cell  = cellfun(Face_Index_Splitter_fun, face_vertex_parts, 'UniformOutput', false);
    
        % Vertex indices
        Vertex_Index_fun = @(vertex_indices) Column_Deal(vertex_indices);             
        
        if Textures == true
            % Texture indices are also given
            [face_vertex_index_matrix, face_texture_index_matrix, face_normal_index_matrix] = cellfun(Vertex_Index_fun, face_vertex_index_cell);
    
        else
            % Otherwise only the vertices and normals
            [face_vertex_index_matrix, face_normal_index_matrix]    = cellfun(Vertex_Index_fun, face_vertex_index_cell);
            face_texture_index_matrix                               = [];
        end
    
        %--% Face centres %--%
        % Divide the indices per face
        [number_faces, number_face_elements]    = size(face_vertex_index_matrix);
        face_vertex_index_cell                  = mat2cell(face_vertex_index_matrix, ones(1, number_faces), number_face_elements);
    
        % The centre is defined as the average of the vertex points
        Face_Centre_fun     = @(face_vertex_index_list) mean(vertex_coord_matrix(face_vertex_index_list, :), 1);
        face_centre_cell    = cellfun(Face_Centre_fun, face_vertex_index_cell, 'UniformOutput', false);
        face_centre_matrix  = vertcat(face_centre_cell{:});
    
        %--% Face normal vectors %--%
        % Computed only if the faces are triangular
        if number_face_elements == 3
            % Function to compute triangle normal vectors
            Face_Normal_fun             = @(first_index, second_index, third_index) cross(vertex_coord_matrix(second_index, :) - vertex_coord_matrix(first_index, :), vertex_coord_matrix(third_index, :) - vertex_coord_matrix(first_index, :));
            face_normal_vector_cell     = arrayfun(Face_Normal_fun, face_vertex_index_matrix(:, 1), face_vertex_index_matrix(:, 2), face_vertex_index_matrix(:, 3), 'UniformOutput', false);
            face_normal_vector_matrix   = vertcat(face_normal_vector_cell{:});
    
            % Ensure they are unit length
            vector_length_list          = sqrt(sum(face_normal_vector_matrix.^2, 2));
            face_normal_vector_matrix   = face_normal_vector_matrix ./ vector_length_list;
        else
            face_normal_vector_matrix = [];
        end
    
        %--% Face structure %--%
        Faces = struct('vertex_indices', face_vertex_index_matrix, 'texture_indices', face_texture_index_matrix, 'normal_indices', face_normal_index_matrix, 'normal_vector', face_normal_vector_matrix, 'centre', face_centre_matrix);

    %% Plot %%
        if Plot == true
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])  
            
            hold on
            grid on
        
            % Vertices
            scatter3(vertex_coord_matrix(:, 1), vertex_coord_matrix(:, 2), vertex_coord_matrix(:, 3), 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', 'DisplayName', 'Vertices');
            
            % Triangular mesh properties
            if number_face_elements == 3
                % Mesh
                trimesh(face_vertex_index_matrix, vertex_coord_matrix(:, 1), vertex_coord_matrix(:, 2), vertex_coord_matrix(:, 3), 'LineWidth', 1, 'FaceColor', 'none', 'EdgeColor', 'k', 'DisplayName', 'Mesh');
        
                % Vector scale
                vertex_UB       = max(vertex_coord_matrix, [], 1);
                vertex_LB       = min(vertex_coord_matrix, [], 1);
                vector_scale    = norm(vertex_UB - vertex_LB) / 10;
        
                % Face normal vectors
                for f = 1 : number_faces
                    face_centre = face_centre_matrix(f, :);
                    face_normal = vector_scale * face_normal_vector_matrix(f, :);           % Scaled
        
                    pl_normal = plot3(face_centre(1) + [0, face_normal(1)], face_centre(2) + [0, face_normal(2)], face_centre(3) + [0, face_normal(3)], 'LineWidth', 2, 'Color', 'c', 'DisplayName', 'Face normal');
        
                    if f > 1
                        pl_normal.HandleVisibility = 'Off';
                    end
                end
            end
        
            % Axes
            xlabel('x [m]');
            ylabel('y [m]');
            zlabel('z [m]');
        
            axis equal
            view(45, 45);
        
            % Legend
            legend('show', 'location', 'eastoutside');
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            % Pause message
            disp('The object file has been read. The figure will close and script end upon a key-press.');
            pause();

            close(1);
        end
end