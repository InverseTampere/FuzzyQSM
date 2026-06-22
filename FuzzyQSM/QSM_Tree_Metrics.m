% This script determines the tree-level metrics given cylinder and branch-order metrics

function Tree_Metrics = QSM_Tree_Metrics(Cylinder_Metrics, Branch_Order_Metrics)

    %% Structure inputs %%
        % Cylinder metrics
        cylinder_axis_matrix        = Cylinder_Metrics.axis;
        cylinder_start_matrix       = Cylinder_Metrics.start;
        cylinder_radius_list        = Cylinder_Metrics.radius;
        cylinder_length_list        = Cylinder_Metrics.length;

        if isfield(Cylinder_Metrics, 'UnmodRadius')
            cylinder_UnmodRadius_list = Cylinder_Metrics.UnmodRadius;
        else
            cylinder_UnmodRadius_list = cylinder_radius_list;
        end

        % Branch-order metrics
        branch_order_volume_list    = Branch_Order_Metrics.volume;
        branch_order_length_list    = Branch_Order_Metrics.length;

        if isfield(Branch_Order_Metrics, 'UnmodVolume')
            branch_order_UnmodVolume_list = Branch_Order_Metrics.UnmodVolume;
        else
            branch_order_UnmodVolume_list = branch_order_volume_list;
        end

    %% Manual inputs %%
        Diagnostics     = false;        % [true, false] Shows the fitted cylinders

    %% Tree metrics %%
        %--% Volume %--%
        % Modified
        total_volume    = sum(branch_order_volume_list);
        stem_volume     = branch_order_volume_list(1);                          % The branch orders are given in order, thus the stem is the first entry
        branch_volume   = sum(branch_order_volume_list(2 : end));

        % Unmodified
        total_UnmodVolume   = sum(branch_order_UnmodVolume_list);
        stem_UnmodVolume    = branch_order_UnmodVolume_list(1);                 % The branch orders are given in order, thus the stem is the first entry
        branch_UnmodVolume  = sum(branch_order_UnmodVolume_list(2 : end));

        %--% DBH %--%
        height_start_list   = cylinder_start_matrix(:, 3);
        height_end_list     = height_start_list + cylinder_length_list .* cylinder_axis_matrix(:, 3);

        height_matrix       = [height_start_list, height_end_list];
        height_start_list   = min(height_matrix, [], 2);                    % In case the cylinder axis points down
        height_end_list     = max(height_matrix, [], 2);

        min_height          = min(height_start_list);
        height_start_list   = height_start_list - min_height;           % Ensure that the tree starts at height 0
        height_end_list     = height_end_list - min_height;

        BH                  = 1.3;
        BH_cylinder_bool    = height_start_list < BH & height_end_list > BH;

        DBH         = 2 * cylinder_radius_list(BH_cylinder_bool);       
        UnmodDBH    = 2 * cylinder_UnmodRadius_list(BH_cylinder_bool);

        if sum(BH_cylinder_bool) > 1
            % The largest diameter is taken and a warning message is displayed
            DBH         = max(DBH);
            UnmodDBH    = max(UnmodDBH);
            warning('%i cylinders were found at breast-height. The largest diameter of %.3g m is used.', sum(BH_cylinder_bool), DBH);
        elseif sum(BH_cylinder_bool) == 0
            % It can happen that there is a slight gap at BH. Then the average is taken of the cylinders below and above it
            cylinder_below = find(height_end_list < BH, 1, 'last');
            cylinder_above = find(height_start_list > BH, 1, 'first');

            DBH         = cylinder_radius_list(cylinder_below) + cylinder_radius_list(cylinder_above);              % Note that the diameter means it's just the sum
            UnmodDBH    = cylinder_UnmodRadius_list(cylinder_below) + cylinder_UnmodRadius_list(cylinder_above);

            % A warning message shows the gap
            gap_size = height_start_list(cylinder_above) - height_end_list(cylinder_below);
            warning('A gap of %.3g m was found between the two cylinders at breast-height.', gap_size);
        end

        %--% Tree height %--%
        cylinder_end_matrix     = cylinder_start_matrix + cylinder_length_list .* cylinder_axis_matrix;
        num_dim                 = size(cylinder_start_matrix, 2);
        cylinder_height_list    = [cylinder_start_matrix(:, num_dim); cylinder_end_matrix(:, num_dim)];
        tree_height             = max(cylinder_height_list);

        %--% Total branch length %--%
        total_branch_length = sum(branch_order_length_list(2 : end));       % The stem is excluded

        %--% Structure %--%
        Tree_Metrics = struct('TotalVolume', total_volume, 'TrunkVolume', stem_volume, 'BranchVolume', branch_volume, 'TotalUnmodVolume', total_UnmodVolume, 'TrunkUnmodVolume', stem_UnmodVolume, 'BranchUnmodVolume', branch_UnmodVolume, 'DBHqsm', DBH, 'UnmodDBHqsm', UnmodDBH, 'TreeHeight', tree_height, 'BranchLength', total_branch_length);

    %% Diagnostics plot %%
        if Diagnostics == true
            % Diagnostics plot
            number_cylinders = length(cylinder_radius_list);
            
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])   

            hold on
            grid on

            for c = 1 : number_cylinders
                cyl_radius  = cylinder_radius_list(c);
                cyl_length  = cylinder_length_list(c);
                cyl_axis    = cylinder_axis_matrix(c, :);
                cyl_start   = cylinder_start_matrix(c, :);
                cyl_centre  = cyl_start + cyl_length/2*cyl_axis;

                if BH_cylinder_bool(c) == true
                    cyl_colour = 'b';
                else
                    cyl_colour = 'r';
                end

                % The cylinder surface
                number_coord = 1e2;
                [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cyl_radius, cyl_length, cyl_centre, cyl_axis, number_coord);
                surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', cyl_colour, 'FaceAlpha', 0.10, 'LineWidth', 2);
            end

            % Axes
            xlabel('x [m]')
            ylabel('y [m]')
            zlabel('z [m]')

            axis equal
            view(45, 45);

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
    
            hold off

            % Pause message
            disp('The tree metrics have been determined. The script ends and figure closes upon a key-press.');
            pause();

            close(1);
        end
end