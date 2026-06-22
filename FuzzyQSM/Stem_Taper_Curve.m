% The circle radii in the tree model's stem are compared to the cylinder radii of the given QSM
% Note that Object_Circle_Geometry's cylinder indices are expected to refer to the given QSM's
% The error is defined as QSM - object

function [height_list, object_radius_list, QSM_radius_list, radius_error_list, radius_rel_error_list] = Stem_Taper_Curve(QSM, Object_Circle_Geometry, Parallel_Pool, Plot)

    %% Inputs %%
        % Cylinder geometry
        cylinder_BO_list        = QSM.cylinder.BranchOrder;

        % Object circles
        obj_circle_BO_list      = Object_Circle_Geometry.BranchOrder;

    %% Stem circle to cylinder association %%
        % Removing non-stem cylinders and circles
        stem_cylinder_bool      = cylinder_BO_list == 0;
        QSM_Stem_Cylinders      = Structure_Boolean(QSM.cylinder, stem_cylinder_bool);
        number_cylinders        = sum(stem_cylinder_bool);
        cylinder_radius_list    = QSM_Stem_Cylinders.radius;

        stem_circle_bool        = obj_circle_BO_list == 0;
        Stem_Circle_Geometry    = Structure_Boolean(Object_Circle_Geometry, stem_circle_bool);
        
        Stem_Circle_Geometry.number_circles = sum(stem_circle_bool);

        % Association
        Stem_QSM.cylinder           = QSM_Stem_Cylinders;
        Stem_Circle_Geometry        = Smoothed_Model_Circle_Cylinder_Association(Stem_QSM, Stem_Circle_Geometry, Parallel_Pool);
        stem_circle_cyl_index_list  = Stem_Circle_Geometry.cylinder_index;
        stem_circle_radius_list     = Stem_Circle_Geometry.radius;
        stem_circle_centre_matrix   = Stem_Circle_Geometry.centre;
        stem_circle_height_list     = stem_circle_centre_matrix(:, 3);

    %% Radii at different heights %%  
        % Data for the stem cylinders
        cylinder_radius_cell    = cell(1, number_cylinders);
        stem_circle_height_cell = cell(1, number_cylinders);
        stem_circle_radius_cell = cell(1, number_cylinders);

        for c = 1 : number_cylinders
            % This cylinder's data
            cylinder_radius = cylinder_radius_list(c);

            % Circles associated with it
            circle_bool                 = stem_circle_cyl_index_list == c;
            stem_circle_height_cell{c}  = stem_circle_height_list(circle_bool);
            stem_circle_radius_cell{c}  = stem_circle_radius_list(circle_bool);

            % The cylinder radius is entered once for each circle
            number_circles          = sum(circle_bool);
            cylinder_radius_cell{c} = cylinder_radius * ones(number_circles, 1);
        end

        % Total data
        height_list             = vertcat(stem_circle_height_cell{:});
        [height_list, order]    = sort(height_list);

        object_radius_list = vertcat(stem_circle_radius_cell{:});
        object_radius_list = object_radius_list(order);

        QSM_radius_list = vertcat(cylinder_radius_cell{:});
        QSM_radius_list = QSM_radius_list(order);

        % Errors
        radius_error_list       = QSM_radius_list - object_radius_list;
        radius_rel_error_list   = radius_error_list ./ object_radius_list * 100;

    %% Plot %%
        if Plot == true
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85]);
            set(gcf, 'color', [1, 1, 1]);

            % First subplot showing both data sets
            subplot(1, 3, 1)
            hold on
            grid on

            plot(object_radius_list, height_list, 'LineWidth', 2, 'color', 'r', 'DisplayName', 'Object');
            plot(QSM_radius_list, height_list, 'LineWidth', 2, 'color', 'b', 'DisplayName', 'QSM');

            xlabel('radius [m]');
            ylabel('height [m]');

            legend('show', 'location', 'northeast');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
    
            hold off

            % Second subplot for the error
            subplot(1, 3, 2)
            hold on
            grid on

            plot(radius_error_list, height_list, 'LineWidth', 2, 'color', 'm');

            xlabel('radius error [m]');
            ylabel('height [m]');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
    
            hold off

            % Third subplot for the relative error
            subplot(1, 3, 3)
            hold on
            grid on

            plot(radius_rel_error_list, height_list, 'LineWidth', 2, 'color', 'c');

            xlabel('radius rel. error [%]');
            ylabel('height [m]');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
    
            hold off

            % Pause message
            disp('The stem taper curve has been determined. The figure will close and script end upon a key-press.');
            pause();

            close(1);
        end
end