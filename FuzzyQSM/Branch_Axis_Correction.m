% If the axis directions between subsequent cylinders deviate more than the given threshold and if it does not for the alternative cylinder fit,
% that alternative axis direction is chosen
% In addition it ensures that the axis directions point down the branch

function Cylinders = Branch_Axis_Correction(Cylinders, Alternative_Cylinders)

    %% Structure inputs %%
        % Cylinder properties
        cyl_axis_matrix         = Cylinders.axis;
        cyl_centre_matrix       = Cylinders.centre;
        cyl_radius_list         = Cylinders.radius;
        cyl_length_list         = Cylinders.length;
        cyl_branch_index_list   = Cylinders.branch;
        cyl_branch_pos_list     = Cylinders.PositionInBranch;

        % Alternative cylinder geometry
        alt_cyl_axis_matrix     = Alternative_Cylinders.axis;

    %% Manual inputs %%
        angle_threshold         = 45;               % [deg] Allowed deviation between cylinder axes
                                                    %       If >90, axes are not substituted but their direction is still ensured to be down the branch
        Diagnostics             = false;            % [true, false] Shows the result for a specific branch

    %% Branch axis alignment %%        
        % Angle threshold in radians
        angle_threshold = deg2rad(angle_threshold);

        % New axes
        [number_cylinders, ~]   = size(cyl_axis_matrix);
        new_cyl_axis_matrix     = cyl_axis_matrix;
        axis_changed_bool       = false(number_cylinders, 1);

        % Unique branches
        branch_index_list   = unique(cyl_branch_index_list);
        number_branches     = length(branch_index_list);

        for b = 1 : number_branches
            % Cylinders belonging to this branch
            branch_index            = branch_index_list(b);
            branch_cylinders        = find(cyl_branch_index_list == branch_index);
            number_branch_cylinders = length(branch_cylinders);

            if number_branch_cylinders == 1        
                continue
            end

            % Ensuring they are placed in the correct order
            branch_cyl_pos_list = cyl_branch_pos_list(branch_cylinders);
            [~, order]          = sort(branch_cyl_pos_list);
            branch_cylinders    = branch_cylinders(order);

            % Their geometry
            branch_cyl_axis_matrix      = cyl_axis_matrix(branch_cylinders, :);
            branch_cyl_centre_matrix    = cyl_centre_matrix(branch_cylinders, :);
            branch_cyl_radius_list      = cyl_radius_list(branch_cylinders);
            branch_cyl_length_list      = cyl_length_list(branch_cylinders);
            alt_branch_cyl_axis_matrix  = alt_cyl_axis_matrix(branch_cylinders, :);

            %--% Ensuring they point down the branch %--%
            % Ensure that the first cylinder axis points towards the next cylinder's centre, as otherwise they can all point up the branch
            first_cyl_centre    = branch_cyl_centre_matrix(1, :);
            first_cyl_axis      = branch_cyl_axis_matrix(1, :);
            second_cyl_centre   = branch_cyl_centre_matrix(2, :);
            [~, delta, ~]       = Point_to_Vector_Projection(second_cyl_centre, first_cyl_axis, first_cyl_centre);

            branch_cyl_axis_matrix(1, :) = sign(delta) * first_cyl_axis;

            % The dot product between subsequent cylinder axes should be positive
            branch_cyl_axis_matrix_prev = branch_cyl_axis_matrix(1 : number_branch_cylinders - 1, :);
            branch_cyl_axis_matrix_next = branch_cyl_axis_matrix(2 : number_branch_cylinders, :);
            subseq_dot_product_list     = dot(branch_cyl_axis_matrix_prev, branch_cyl_axis_matrix_next, 2);

            branch_cyl_axis_matrix(2 : number_branch_cylinders, :) = sign(subseq_dot_product_list) .* branch_cyl_axis_matrix(2 : number_branch_cylinders, :);

            % The dot products between these and the alternate axes should also be positive
            alt_dot_product_list        = dot(branch_cyl_axis_matrix, alt_branch_cyl_axis_matrix, 2);
            alt_branch_cyl_axis_matrix  = sign(alt_dot_product_list) .* alt_branch_cyl_axis_matrix;

            %--% Checking alignment in order %--%
            new_branch_cyl_axis_matrix  = branch_cyl_axis_matrix;
            changed_bool                = false(number_branch_cylinders, 1); 

            % Intermediate cylinders
            for c = 2 : number_branch_cylinders - 1       
                % Alignment between current and the previous and next cylinders
                cyl_axis_prev   = branch_cyl_axis_matrix(c - 1, :);
                cyl_axis        = branch_cyl_axis_matrix(c, :);
                cyl_axis_next   = branch_cyl_axis_matrix(c + 1, :);
                
                [~, ~, aligned_prev, ~]  = Vector_Alignment_Check(cyl_axis, cyl_axis_prev, angle_threshold);
                [~, ~, aligned_next, ~]  = Vector_Alignment_Check(cyl_axis, cyl_axis_next, angle_threshold);

                % If neither is aligned, the alternative axis is checked
                if ~aligned_prev && ~aligned_next
                    % Check if the alternative axis is aligned
                    alt_cyl_axis            = alt_branch_cyl_axis_matrix(c, :);
                    [~, ~, aligned_prev, ~] = Vector_Alignment_Check(alt_cyl_axis, cyl_axis_prev, angle_threshold);
                    [~, ~, aligned_next, ~] = Vector_Alignment_Check(alt_cyl_axis, cyl_axis_next, angle_threshold);

                    if aligned_prev && aligned_next
                        new_branch_cyl_axis_matrix(c, :)    = alt_cyl_axis;
                        changed_bool(c)                     = true;
                    end
                end
            end

            % First cylinder can only be compared to the second
            first_cyl_axis  = branch_cyl_axis_matrix(1, :);
            second_cyl_axis = branch_cyl_axis_matrix(2, :);

            [~, ~, aligned_first, ~] = Vector_Alignment_Check(first_cyl_axis, second_cyl_axis, angle_threshold);

            if ~aligned_first
                first_alt_cyl_axis          = alt_branch_cyl_axis_matrix(1, :);
                [~, ~, aligned_first, ~]    = Vector_Alignment_Check(first_alt_cyl_axis, second_cyl_axis, angle_threshold);

                if aligned_first
                    new_branch_cyl_axis_matrix(1, :) = first_alt_cyl_axis;
                end
            end

            % Last cylinder can only be compared to the second-to-last
            last_cyl_axis           = branch_cyl_axis_matrix(number_branch_cylinders, :);
            second_to_last_cyl_axis = branch_cyl_axis_matrix(number_branch_cylinders - 1, :);

            [~, ~, aligned_last, ~] = Vector_Alignment_Check(last_cyl_axis, second_to_last_cyl_axis, angle_threshold);

            if ~aligned_last
                last_alt_cyl_axis       = alt_branch_cyl_axis_matrix(number_branch_cylinders, :);
                [~, ~, aligned_last, ~] = Vector_Alignment_Check(last_alt_cyl_axis, second_to_last_cyl_axis, angle_threshold);

                if aligned_last
                    new_branch_cyl_axis_matrix(number_branch_cylinders, :) = last_alt_cyl_axis;
                end
            end

            % Update the complete arrays
            new_cyl_axis_matrix(branch_cylinders, :)    = new_branch_cyl_axis_matrix;
            axis_changed_bool(branch_cylinders)         = changed_bool;

            %--% Diagnostics plot %--%
            if Diagnostics == true
                %--% Plot %--%
                % Plotting inputs
                cyl_cmap        = cbrewer('qual', 'Set1', number_branch_cylinders);         % Colours
                cyl_cmap        = max(cyl_cmap, 0);
                cyl_cmap        = min(cyl_cmap, 1);
                number_coord    = 1e3;                                                      % Number of points per cylinder
            
                figure(1)
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])     
            
                % First subplot for the old cylinders
                subplot(1, 2, 1)
                title('Old cylinders');
                hold on
                grid on

                % Second for the new cylinders
                subplot(1, 2, 2)
                title('New cylinders');
                hold on
                grid on
            
                % Cylinders in order
                for c = 1 : number_branch_cylinders
                    cylinder_colour = cyl_cmap(c, :);
            
                    % Old cylinders
                    cylinder_centre = branch_cyl_centre_matrix(c, :);
                    cylinder_radius = branch_cyl_radius_list(c);
                    cylinder_axis   = branch_cyl_axis_matrix(c, :);
                    cylinder_length = branch_cyl_length_list(c);
            
                    subplot(1, 2, 1)
                    [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, cylinder_length, cylinder_centre, cylinder_axis, number_coord);
                    surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', cylinder_colour, 'FaceAlpha', 0.25, 'LineWidth', 2);
            
                    % New cylinder
                    new_cylinder_axis = new_branch_cyl_axis_matrix(c, :);
            
                    subplot(1, 2, 2)
                    [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, cylinder_length, cylinder_centre, new_cylinder_axis, number_coord);
                    surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', cylinder_colour, 'FaceAlpha', 0.25, 'LineWidth', 2);        
                end
            
                % Axes
                subplot(1, 2, 1)
                xlabel('x [m]')
                ylabel('y [m]')
                zlabel('z [m]')
            
                axis equal
                view(45, 45);
                                
                subplot(1, 2, 2)
                xlabel('x [m]')
                ylabel('y [m]')
                zlabel('z [m]')
            
                axis equal
                view(45, 45);
                        
                % Font size
                subplot(1, 2, 1)
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
            
                hold off

                subplot(1, 2, 2)
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
            
                hold off

                % Pause message
                fprintf('%i axes have been corrected. Upon a key-press the figure closes and the script continues. \n', sum(changed_bool));
                pause();
                
                close(1);
            end
        end

        % Updating the structure
        Cylinders.axis          = new_cyl_axis_matrix;
        Cylinders.AxisChanged   = axis_changed_bool;

end