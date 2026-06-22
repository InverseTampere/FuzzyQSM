% The given point is projected onto the ellipse using the shortest path

function [projected_point, distance] = Point_to_Ellipse_Projection(ellipse_centre, ellipse_radii, ellipse_axes, point, Plot, Diagnostics)

    %% Inputs %%
        ellipticity_threshold   = 1e-4;     % [-] For the ellipse to be considered a circle
        numerical_margin        = 1e-6;     % [-] Margin taken into account due to rounding errors

    %% Ellipticity check %%
        % The scheme can break if the radii are highly equal as it revolves around solving a quartic which is then in truth a quadratic polynomial of a different form
        ellipticity             = abs(ellipse_radii(1)/ellipse_radii(2) - 1);           % If the radii are equal, this reduces to zeros
    
        % If it is not elliptic, the point is simply projected onto the ellipse as if it were a circle
        if ellipticity < ellipticity_threshold
            % The projected point is placed on the radius
            point_t         = point - ellipse_centre;
            proj_point_t    = mean(ellipse_radii) * point_t / norm(point_t);       
            projected_point = proj_point_t + ellipse_centre;
    
        % Otherwise it is treated as an ellipse and solved accordingly
        else

    %% Coordinate transformation %%
        % An ellipse-centered, ellipse-oriented coordinate frame is used
        point_t = point - ellipse_centre;
        point_r = (ellipse_axes * point_t')';
        
        % The point is placed in the first quarter
        sign_list   = sign(point_r);
        point_q     = sign_list .* point_r;
    
    %% Projection %%
        % If the point is at the origin, the solution in the first quarter equals the positive minor axis
        num_dim = length(point_q);
        abs_q   = abs(point_q);
        
        if max(abs_q) < numerical_margin            
            [minor_radius, minor_axis_bool] = min(ellipse_radii);
            projected_point_q               = minor_radius * minor_axis_bool;
            
            sign_list = ones(1, num_dim);       % Note that the sign of zero equals zero, and so for the point to be retransformed the signs are changed to 1

        % If one of the coordinates is zero, the projected point is found more easily
        elseif min(abs_q) < numerical_margin && max(abs_q) > numerical_margin
            % The non-zero dimension
            nonzero_bool = abs_q > numerical_margin;

            % Projected point
            t = -ellipse_radii(nonzero_bool)^2 + ellipse_radii(nonzero_bool)*point_q(nonzero_bool);
            projected_point_q = point_q ./ (1 + t ./ ellipse_radii.^2);

        else
            % Coefficients
            a = -1;
            b = -2*sum(ellipse_radii.^2);
            c = sum(ellipse_radii.^2 .* point_q.^2) - sum(ellipse_radii.^4) - 4*prod(ellipse_radii.^2);
            d = 2*prod(ellipse_radii.^2)*(sum(point_q.^2) - sum(ellipse_radii.^2));
            e = prod(ellipse_radii.^2)*(ellipse_radii(2)^2*point_q(1)^2 + ellipse_radii(1)^2*point_q(2)^2 - prod(ellipse_radii.^2));

            [t_roots, ~, ~] = Quartic_Polynomial(a, b, c, d, e, Diagnostics, Diagnostics, Diagnostics);

            % The solution to t is the largest root
            t = max(t_roots);
            
            % Projected point
            projected_point_q = point_q ./ (1 + t ./ ellipse_radii.^2);
        end

        %% Transformation to original coordinate system %%
            % No longer restricted to the first quarter
            projected_point_r   = sign_list .* projected_point_q;
            
            % Rotated and translated
            projected_point_t   = projected_point_r * ellipse_axes;
            projected_point     = projected_point_t + ellipse_centre; 
        end
    
        % The distance between the point and its projection
        distance = norm(projected_point - point);

    %% Plot %%
        if Plot == true
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])     

            hold on
            grid on

            % Ellipse
            number_coord            = 1e3;
            ellipse_coord_matrix    = Ellipse_Coordinate_Generator(ellipse_centre, ellipse_axes, ellipse_radii, number_coord);    
            plot(ellipse_coord_matrix(:, 1), ellipse_coord_matrix(:, 2), 'LineWidth', 2, 'color', 'r', 'DisplayName', 'Ellipse');

            % Point
            scatter(point(1), point(2), 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'Point');

            % Projected point
            scatter(projected_point(1), projected_point(2), 'filled' ,'MarkerFaceColor', 'm', 'DisplayName', 'Projected point');

            % Axes
            axis equal
            xlabel('x [m]');
            ylabel('y [m]');

            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    

            % Pause
            disp('The point has been projected. The script will finish and figure will close upon a key-press');
            pause();
            
            close(1);            
        end
end