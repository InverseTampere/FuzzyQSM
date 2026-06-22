% The smoothed tree model consists of segments between circles
% Using the mesh between each connected circle, the volume and min/max radius of each segment is determined

function Segment_Metrics = Smoothed_Model_Segment_Metrics(Mesh, Object_Circle_Geometry, Original_QSM, Parallel_Pool)

    %% Structure inputs %%
        % Mesh
        vertex_matrix                   = Mesh.Vertices.coordinates;
        triangle_vertex_index_matrix    = Mesh.Triangles.vertex_indices;

        % Object circle geometry
        circle_vertex_start_index_list  = Object_Circle_Geometry.vertex_start_indices;
        circle_vertex_end_index_list    = Object_Circle_Geometry.vertex_end_indices;
        circle_branch_index_list        = Object_Circle_Geometry.branch_index;
        number_circles                  = Object_Circle_Geometry.number_circles;

        % Original QSM
        branch_order_list               = Original_QSM.branch.order;

        % Parallel pool
        Parallel_Loop                   = Parallel_Pool.Parallel_Loop;
        max_cores                       = Parallel_Pool.max_cores;
        idle_timeout                    = Parallel_Pool.idle_timeout;

    %% Initiate the parallel pool %%
        if Parallel_Loop == true
            % The parallel pool is started
            if isempty(max_cores)
                number_cores = feature('numcores');
            else
                number_cores = max_cores;
            end
            
            Parallel_Pool_Starter(idle_timeout, number_cores)
        else
            % Otherwise the number of cores used in parallel are set to 0 to make the parfor run as a for
            number_cores = 0;
        end

    %% Segment volumes %%
        % Computed one circle at a time
        segment_volume_cell         = cell(1, number_circles);
        segment_circles_index_cell  = cell(1, number_circles);

        DQ      = parallel.pool.DataQueue;
        tick    = 0;
        N       = number_circles;
        afterEach(DQ, @ProgressUpdate);
    
        parfor (i = 1 : number_circles, number_cores)
            % The circle's data
            circle_i_start_index    = circle_vertex_start_index_list(i);
            circle_i_end_index      = circle_vertex_end_index_list(i);
    
            % To close the mesh, the circle's vertices have to be triangulated
            number_vertices_i                     = circle_i_end_index - circle_i_start_index + 1;
            circle_i_triangle_vertex_index_matrix = circle_i_start_index + [zeros(number_vertices_i - 2, 1), (1 : number_vertices_i - 2)', (2 : number_vertices_i - 1)'];
            
            % Already existing triangles which contain the circle's vertices form the flanks
            flank_i_triangle_bool                 = any(triangle_vertex_index_matrix >= circle_i_start_index & triangle_vertex_index_matrix <= circle_i_end_index, 2);
            flank_i_triangle_vertex_index_matrix  = triangle_vertex_index_matrix(flank_i_triangle_bool, :);
            
            % It is connected to any circles which are at least partially within the range of vertex indices
            [min_flank_i_vertex_index, max_flank_i_vertex_index] = deal(min(flank_i_triangle_vertex_index_matrix(:)), max(flank_i_triangle_vertex_index_matrix(:)));
            
            connected_circles = find(circle_vertex_start_index_list >= min_flank_i_vertex_index & circle_vertex_start_index_list <= max_flank_i_vertex_index);          %#ok<PFBNS>
            
            % To avoid duplicate computations, only connect circles past this one are considered
            connected_circles(connected_circles <= i)   = [];
            number_connected_circles                    = length(connected_circles);
            
            % Volumes in the segment between the current circle and each connected circle
            circle_i_segment_volume_list = NaN(number_connected_circles, 1);
    
            for j = 1 : number_connected_circles
                % The connected circle's vertices
                connected_circle        = connected_circles(j);
                circle_j_start_index    = circle_vertex_start_index_list(connected_circle);
                circle_j_end_index      = circle_vertex_end_index_list(connected_circle);       %#ok<PFBNS>
    
                % To close the mesh, the circle's vertices have to be triangulated
                number_vertices_j                       = circle_j_end_index - circle_j_start_index + 1;
                circle_j_triangle_vertex_index_matrix   = circle_j_start_index + [zeros(number_vertices_j - 2, 1), (1 : number_vertices_j - 2)', (2 : number_vertices_j - 1)'];
    
                % The mesh connecting the two
                flank_ij_triangle_bool                  = any(flank_i_triangle_vertex_index_matrix >= circle_j_start_index & flank_i_triangle_vertex_index_matrix <= circle_j_end_index, 2);
                flank_ij_triangle_vertex_index_matrix   = flank_i_triangle_vertex_index_matrix(flank_ij_triangle_bool, :);
    
                % Watertight segment mesh
                segment_ij_triangle_vertex_index_matrix = [circle_i_triangle_vertex_index_matrix; circle_j_triangle_vertex_index_matrix; flank_ij_triangle_vertex_index_matrix];
    
                segment_ij_mesh = surfaceMesh(vertex_matrix, segment_ij_triangle_vertex_index_matrix);            
                watertightness  = isWatertight(segment_ij_mesh);
    
                if ~watertightness                      % The mesh is not necessarily watertight, for instance when the circles intersect or at a bifurcation
                    number_connected_circles = number_connected_circles - 1;
                    continue
                end
    
                % Its volume
                segment_ij_volume               = Triangular_Mesh_Volume(segment_ij_triangle_vertex_index_matrix, vertex_matrix);
                circle_i_segment_volume_list(j) = segment_ij_volume;
            end
    
            % Adding this circle's segments to the cell array
            NaN_bool = isnan(circle_i_segment_volume_list);                         % A watertight mesh could not be constructed

            segment_volume_cell{i}          = circle_i_segment_volume_list(~NaN_bool);
            segment_circles_index_matrix    = [repmat(i, [number_connected_circles, 1]), connected_circles(~NaN_bool)];
            segment_circles_index_cell{i}   = segment_circles_index_matrix;
                
            % Progress update
            send(DQ, i);
        end

        % Combined arrays
        segment_volume_list             = vertcat(segment_volume_cell{:});
        segment_circles_index_matrix    = vertcat(segment_circles_index_cell{:});

        % Branch indices and orders of the segments
        segment_branch_index_list = circle_branch_index_list(segment_circles_index_matrix(:, 1));
        segment_branch_order_list = branch_order_list(segment_branch_index_list);

        % Output structure
        Segment_Metrics = struct('volume', segment_volume_list, 'circle_indices', segment_circles_index_matrix, 'branch_index', segment_branch_index_list, 'BranchOrder', segment_branch_order_list);

    %% Local functions %%
        % Progress update
        function ProgressUpdate(~)
            % Ensures that at most every P percent is printed
            P = 1;
    
            % The tick is updated
            tick = tick + 1;    
    
            % The last tick's and current tick's progress relative to P
            progress_last   = floor((tick - 1) / N * 100 / P);
            progress        = floor(tick / N * 100 / P);
    
            if progress - progress_last >= 1
                fprintf('   Segment volume calculation progress: %i %% \n', progress);
            end            
        end
end