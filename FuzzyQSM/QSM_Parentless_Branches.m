% The QSM can have branches of order 0 which are not part of the stem because they lack a parent
% Here they (and child branches) are given a NaN branch order

function [QSM, disconnected_branches] = QSM_Parentless_Branches(QSM)

    %% Structure inputs %%
        % Branches
        Branches                = QSM.branch;
        branch_order_list       = double(Branches.order);
        branch_parent_list      = double(Branches.parent);

        % Cylinders
        Cylinders               = QSM.cylinder;
        cylinder_branch_list    = double(Cylinders.branch);
        
    %% Manual inputs %%
        Diagnostics             = false;         % [true, false] Shows the child branches being found at every step
        Plot                    = false;         % [true, false] Shows which branches were found to be disconnected at the end

   %% Detecting disconnected branches %%
        % They are order 0 but are not the stem, which is the first branch
        zero_order_bool = branch_order_list == 0;

        % As a check, their parent should be zero
        no_parent_bool  = branch_parent_list == 0;

        % The result is the combined boolean, with the first branch (stem) removed
        disconnected_branches = find(zero_order_bool & no_parent_bool);
        disconnected_branches = setdiff(disconnected_branches, 1);

        % Change their branch orders to NaN
        branch_order_list(disconnected_branches) = NaN;
        
    %% Children of these branches %%
        % Recursively check if any branches have a disconnected branch as parent
        disconnected_branch_parents = disconnected_branches;

        while ~isempty(disconnected_branch_parents)
            % Cylinders which have a disconnected parent as their parent
            Parent_fun                          = @(branch_parent) ismember(branch_parent, disconnected_branch_parents);
            disconnected_branch_children_bool   = arrayfun(Parent_fun, branch_parent_list);
            disconnected_branch_children        = find(disconnected_branch_children_bool);

            % Their branch orders are specified as NaN
            branch_order_list(disconnected_branch_children_bool) = NaN;

            % Diagnostics plot
            if Diagnostics == true
                % Number of coordinates per cylinder
                number_coord = 1e2;
            
                figure(1)
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])    
            
                hold on
                grid on
            
                number_cylinders = length(cylinder_branch_list);
                for c = 1 : number_cylinders
                    % Cylinder data
                    cylinder_branch = cylinder_branch_list(c);
                    cylinder_radius = Cylinders.radius(c);
                    cylinder_length = Cylinders.length(c);
                    cylinder_axis   = Cylinders.axis(c, :);
                    cylinder_start  = Cylinders.start(c, :);
                    cylinder_centre = cylinder_start + cylinder_length/2 * cylinder_axis;
            
                    % Cylinder colour
                    if ismember(cylinder_branch, disconnected_branch_parents)
                        cylinder_colour = 'b';
                        cylinder_label  = 'Disconnected parent';
                        cylinder_alpha  = 0.5;
                    elseif ismember(cylinder_branch, disconnected_branch_children)
                        cylinder_colour = 'r';
                        cylinder_label  = 'Disconnected child';
                        cylinder_alpha  = 0.5;
                    else
                        cylinder_colour = 'k';
                        cylinder_label  = 'Regular cylinder';
                        cylinder_alpha  = 0.25;
                    end
            
                    % Surface plot
                    [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, cylinder_length, cylinder_centre, cylinder_axis, number_coord);
                    surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', cylinder_colour, 'FaceAlpha', cylinder_alpha, 'LineWidth', 2, 'DisplayName', cylinder_label);
                end
            
                % Axes
                xlabel('x [m]');
                ylabel('y [m]');
                zlabel('z [m]');
                axis equal
                view(45, 45);

                % Legend
                Graphical_Objects = findobj(gca, '-property', 'DisplayName');
                display_names = string(get(Graphical_Objects, 'DisplayName'));
                
                [unique_display_names, unique_indices] = unique(display_names, 'stable');
                legend(Graphical_Objects(unique_indices), unique_display_names, 'location', 'eastoutside');
            
                % Formatting
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
                
                hold off

                % Pause message
                disp('The children of these disconnected branches have been determined. The figure will close and script continue upon a key-press.');
                pause();

                close(1);
            end
        
            % For the next round they are considered parents
            disconnected_branch_parents = disconnected_branch_children;
        end

    %% Update the structure %%
        % The list is updated to include NaN for disconnected branches
        disconnected_branches   = find(isnan(branch_order_list));
        QSM.branch.order        = branch_order_list;

    %% Plot %%
        if Plot == true
            % Number of coordinates per cylinder
            number_coord = 1e2;
        
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
        
            hold on
            grid on
        
            number_cylinders = length(cylinder_branch_list);
            for c = 1 : number_cylinders
                % Cylinder data
                cylinder_branch = cylinder_branch_list(c);
                cylinder_radius = Cylinders.radius(c);
                cylinder_length = Cylinders.length(c);
                cylinder_axis   = Cylinders.axis(c, :);
                cylinder_start  = Cylinders.start(c, :);
                cylinder_centre = cylinder_start + cylinder_length/2 * cylinder_axis;
        
                % Cylinder colour
                if ismember(cylinder_branch, disconnected_branches)
                    cylinder_colour = 'b';
                    cylinder_label  = 'Disconnected cylinder';
                    cylinder_alpha  = 0.5;
                else
                    cylinder_colour = 'k';
                    cylinder_label  = 'Regular cylinder';
                    cylinder_alpha  = 0.25;
                end
        
                % Surface plot
                [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, cylinder_length, cylinder_centre, cylinder_axis, number_coord);
                surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', cylinder_colour, 'FaceAlpha', cylinder_alpha, 'LineWidth', 2, 'DisplayName', cylinder_label);
            end
        
            % Axes
            xlabel('x [m]');
            ylabel('y [m]');
            zlabel('z [m]');
            axis equal
            view(45, 45);

            % Legend
            Graphical_Objects = findobj(gca, '-property', 'DisplayName');
            display_names = string(get(Graphical_Objects, 'DisplayName'));
            
            [unique_display_names, unique_indices] = unique(display_names, 'stable');
            legend(Graphical_Objects(unique_indices), unique_display_names, 'location', 'eastoutside');
        
            % Formatting
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
            
            hold off

            % Pause message
            disp('The disconnected branches have been determined. The figure will close and script end upon a key-press.');
            pause();

            close(1);
        end

end