% The length of the cylinders within the branches are fitted dependently, i.e. by taking subsequent cylinders into account
% The interface is placed using the intersection between the parent cylinder's top plane and extension cylinder's axis
% Note that it presumes the cylinder axes to point down the branch

function Cylinders = Branch_Interface_Fitting(Cylinders)

    %% Structure inputs %%
        % Cylinder properties
        cyl_axis_matrix         = Cylinders.axis;
        cyl_radius_list         = Cylinders.radius;
        cyl_centre_matrix       = Cylinders.centre;
        cyl_length_list         = Cylinders.length;
        cyl_branch_index_list   = Cylinders.branch;
        cyl_branch_pos_list     = Cylinders.PositionInBranch;

    %% Manual inputs %%
        angle_threshold         = 80;               % [deg] If the angle between the axis and plane is greater than this, the interface is not updated 
        Diagnostics             = false;            % [true, false] Shows the result for a specific branch

    %% Branch interface fitting %%        
        % Degrees to radians
        angle_threshold = deg2rad(angle_threshold);

        % New centre locations and lengths
        new_cyl_centre_matrix   = cyl_centre_matrix;
        new_cyl_length_list     = cyl_length_list;

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
            branch_cyl_centre_matrix    = cyl_centre_matrix(branch_cylinders, :);
            branch_cyl_axis_matrix      = cyl_axis_matrix(branch_cylinders, :);
            branch_cyl_length_list      = cyl_length_list(branch_cylinders);
            branch_cyl_radius_list      = cyl_radius_list(branch_cylinders);

            % Cylinder top and bottoms
            branch_cyl_top_matrix = branch_cyl_centre_matrix + branch_cyl_length_list/2 .* branch_cyl_axis_matrix;
            branch_cyl_bot_matrix = branch_cyl_centre_matrix - branch_cyl_length_list/2 .* branch_cyl_axis_matrix;

            %--% Interface fitting %--%
            % Define the extension's bottom based on the parent's top
            new_branch_cyl_bot_matrix = branch_cyl_bot_matrix;

            for c = 1 : number_branch_cylinders - 1        
                % Current cylinder's geometry
                cyl_top     = branch_cyl_top_matrix(c, :);
                cyl_axis    = branch_cyl_axis_matrix(c, :);

                % Next cylinder's geometry
                next_cyl_top    = branch_cyl_top_matrix(c + 1, :);
                next_cyl_axis   = branch_cyl_axis_matrix(c + 1, :);

                % Intersection between next cylinder's axis and current cylinder's top plane
                [next_cyl_bot, ~, incidence_angle, ~] = Plane_Vector_Intersection(cyl_top, cyl_axis, next_cyl_top, next_cyl_axis);

                % If the incidence angle is very high, an intersection could not be computed properly
                if incidence_angle < angle_threshold
                    new_branch_cyl_bot_matrix(c + 1, :) = next_cyl_bot;
                end
            end

            % Resulting centres and lengths
            new_branch_cyl_centre_matrix                = (new_branch_cyl_bot_matrix + branch_cyl_top_matrix) / 2; 
            new_cyl_centre_matrix(branch_cylinders, :)  = new_branch_cyl_centre_matrix;

            new_branch_cyl_length_list              = sqrt(sum((new_branch_cyl_bot_matrix - branch_cyl_top_matrix).^2, 2));
            new_cyl_length_list(branch_cylinders)   = new_branch_cyl_length_list;

            %--% Diagnostics plot %--%
            if Diagnostics == true
                %--% Table %--%
                length_matrix   = [branch_cyl_length_list, new_branch_cyl_length_list];
                number_digits   = 4;
                output_format   = 'Float';
                row_names       = string(1 : number_branch_cylinders);
                column_names    = {'L old [m]', 'L new [m]'};
                Print           = true;
                
                Table_Formatter(length_matrix, number_digits, output_format, row_names, column_names, [], Print);

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
            
                    [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, cylinder_length, cylinder_centre, cylinder_axis, number_coord);
                    old_cyl_surf = surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Old cylinders');
            
                    % New cylinder
                    new_cyl_length = new_cyl_length_list(c);
                    new_cyl_centre = new_cyl_centre_matrix(c, :);
            
                    [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, new_cyl_length, new_cyl_centre, cylinder_axis, number_coord);
                    surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', cylinder_colour, 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', sprintf('New cylinder %i', c));        
            
                    if c > 1
                        old_cyl_surf.HandleVisibility = 'Off';
                    end
                end
            
                % Axes
                xlabel('x [m]')
                ylabel('y [m]')
                zlabel('z [m]')
            
                axis equal
                view(45, 45);
            
                % Legend
                legend('show', 'location', 'eastoutside');
            
                % Font size
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
            
                hold off

                % Pause message
                disp('The branch interfaces have been fitted. Upon a key-press the figure closes and the script continues.');
                pause();
                
                close(1);
            end
        end

        % Updating the structure
        new_cyl_start_matrix    = new_cyl_centre_matrix - new_cyl_length_list/2 .* cyl_axis_matrix;
        new_cyl_volume_list     = pi*cyl_radius_list.^2 .* new_cyl_length_list;

        Cylinders.centre    = new_cyl_centre_matrix;
        Cylinders.length    = new_cyl_length_list;
        Cylinders.volume    = new_cyl_volume_list;
        Cylinders.start     = new_cyl_start_matrix;
end