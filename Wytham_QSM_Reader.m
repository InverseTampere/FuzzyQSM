% This script reads the Wyntham QSM data and creates a QSM structure for it in the format of TreeQSM 2.4

function QSM = Wytham_QSM_Reader(QSM_file_loc, Plot)
    
    %% Create the QSM structure %%
        % QSM data in the old format
        Wyntham_QSM         = load(QSM_file_loc);
        cyl_start_matrix    = Wyntham_QSM.Sta;
        cyl_axis_matrix     = Wyntham_QSM.Axe;
        cyl_length_list     = Wyntham_QSM.Len;
        cyl_radius_list     = Wyntham_QSM.Rad;

        cyl_branch_index_list   = Wyntham_QSM.BoC(:, 1);
        cyl_branch_order_list   = Wyntham_QSM.BoC(:, 2);
        cyl_branch_pos_list     = Wyntham_QSM.BoC(:, 3);
        
        branch_length_list      = Wyntham_QSM.BLen;
        branch_order_list       = Wyntham_QSM.BOrd;
        branch_volume_list      = Wyntham_QSM.BVol;

        % The lowest cylinder should start at the origin
        height_list     = cyl_start_matrix(:, 3);
        [~, bottom_ind] = min(height_list);

        tree_bottom         = cyl_start_matrix(bottom_ind, :);
        cyl_start_matrix_c  = cyl_start_matrix - tree_bottom;
 
        % Cylinder structure
        QSM_Cylinders   = struct('start', cyl_start_matrix_c, 'axis', cyl_axis_matrix, 'length', cyl_length_list, 'radius', cyl_radius_list, 'branch', cyl_branch_index_list, 'BranchOrder', cyl_branch_order_list, 'PositionInBranch', cyl_branch_pos_list);

        % Branch structure
        QSM_Branches    = struct('length', branch_length_list, 'order', branch_order_list, 'volume', branch_volume_list);

        % Tree data
        stem_bool       = branch_order_list == 0;
        stem_volume     = sum(branch_volume_list(stem_bool));
        branch_volume   = sum(branch_volume_list(~stem_bool));
        total_volume    = stem_volume + branch_volume;

        height_start_list   = cyl_start_matrix_c(:, 3);
        height_end_list     = height_start_list + cyl_length_list.* cyl_axis_matrix(:, 3);

        BH                  = 1.3;
        BH_cylinder_bool    = height_start_list < BH & height_end_list > BH;

        if sum(BH_cylinder_bool) ~= 1
            % Plot the QSM
            QSM_Plot();

            % Error message
            error('%i cylinders were found at breast-height.', sum(BH_cylinder_bool));
        end

        DBH = 2 * cyl_radius_list(BH_cylinder_bool);     

        QSM_Tree = struct('TotalVolume', total_volume, 'BranchVolume', branch_volume, 'TrunkVolume', stem_volume, 'DBHqsm', DBH);

        % Full structure
        QSM = struct('cylinder', QSM_Cylinders, 'branch', QSM_Branches, 'treedata', QSM_Tree);

    %% Plot %%
        if Plot == true
            % Plot the QSM
            QSM_Plot();

            % Pause message
            disp('The QSM has been retrieved. Upon a key-press the figure closes and script ends.');
            pause();

            close(1);
        end

        function QSM_Plot()
            % QSM plot inputs
            figure_number   = 1;
            alpha_value     = 1.0;
            number_facets   = 10;

            % Figure
            Fig = figure(figure_number);  
            set(Fig, 'name', QSM_file_loc, 'NumberTitle', 'off');

            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    

            % QSM
            plot_cylinder_model(QSM_Cylinders, "order", figure_number, number_facets, alpha_value);

            % Model
            point_cloud_matrix      = Wyntham_QSM.P;
            point_cloud_matrix_c    = point_cloud_matrix - tree_bottom;

            scatter3(point_cloud_matrix_c(:, 1), point_cloud_matrix_c(:, 2), point_cloud_matrix_c(:, 3), 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');

            % Axes
            xlabel('x [m]');
            ylabel('y [m]');
            zlabel('z [m]');

            view(45, 45);

            axis equal

            % Formatting
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off
        end
end