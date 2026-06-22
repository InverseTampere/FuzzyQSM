% This script fits a circle to the given 2D data by using the least-squares approach
% It requires x and y coordinates to be given, and spits out the centre of the circle (circle_x, circle_y) as well as its radius (circle_r)

% Three methods can be used: algebraic least-squares (Algebraic), geometric least-squares (Geometric) and Pratt's least-squares (Pratt)
% They can be weighted by uncertainty if Weighting is true. This requires the locations of the scanners, laser beam divergence and range uncertainty of the scanner itself
% Alternatively, they can be simply weighted by distance too
% Note however that weighting is not currently implemented for geometric least-squares

function [circle_centre, circle_radius] = Least_Squares_Circle_Fitting(point_cloud_matrix, weight_list, Method, Plot)
    
    %% Inputs %%
        % Threshold for convergence
        min_step_threshold  = 1e-6;         % [%] If the step size is below this value, it has converged
        max_iterations      = 1e2;          % [-] If the maximum number of iterations is exceeded the last values will be used and a message is printed
        
    %% Point cloud dimensions %%
        % For clarity the x and y coordinates of the point cloud have separate variables
        [x_list, y_list]    = Column_Deal(point_cloud_matrix);
        number_points       = length(x_list);

    %% Algebraic least squares %%
        if strcmp(Method, 'Algebraic')
            % The least squares fit matrices (A*x = B)
            A = [x_list, y_list, ones(number_points, 1)];
            B = -(x_list.^2 + y_list.^2);        
            
            % Weighted
            W = [weight_list, weight_list, weight_list];
    
            C = A' * (W .* A);
            D = A' * (W .* B);
    
            % The least squares solution
            circle_parameters_w = C \ D;
    
            circle_x = -circle_parameters_w(1) / 2;
            circle_y = -circle_parameters_w(2) / 2;

            circle_centre = [circle_x, circle_y];
            circle_radius = sqrt(circle_x^2 + circle_y^2 - circle_parameters_w(3));
        
    %% Geometrical least squares %%
        elseif strcmp(Method, 'Geometric')
            %--% An unweighted fit is performed for the initial guess %--%
            % The least squares fit matrices (A*x = B)
            A = [x_list, y_list, ones(number_points, 1)];
            B = -(x_list.^2 + y_list.^2);        
            
            % The least squares solution
            circle_parameters = A \ B;
    
            circle_x        = -circle_parameters(1) / 2;
            circle_y        = -circle_parameters(2) / 2;
            circle_radius   = sqrt(circle_x^2 + circle_y^2 - circle_parameters(3));        
            
            % Gauss-Newton method to solve the non-linear least-squares problem
            convergence = false;
            iter        = 0;
            
            while convergence == false && iter < max_iterations
                iter = iter + 1;
                
                % The residuals
                R = sqrt((x_list - circle_x).^2 + (y_list - circle_y).^2) - circle_radius;
    
                % The derivatives
                dRdx_g = (circle_x - x_list) ./ sqrt((x_list - circle_x).^2 + (y_list - circle_y).^2);
                dRdy_g = (circle_y - y_list) ./ sqrt((x_list - circle_x).^2 + (y_list - circle_y).^2);
                dRdr_g = -ones(number_points, 1);
    
                % The Jacobian
                J = [dRdx_g, dRdy_g, dRdr_g];
    
                % The next geometry parameters
                Beta        = [circle_x; circle_y; circle_radius];
                Beta_new    = Beta - (J' * J)\(J' * R);
    
                circle_x        = Beta_new(1);
                circle_y        = Beta_new(2);
                circle_centre   = [circle_x, circle_y];
                circle_radius   = Beta_new(3);
    
                % Check for convergence of the squared sum of residuals
                sum_residuals = sum(R.^2);
                
                if iter == 1
                    sum_residuals_init = sum_residuals;
                end
    
                R_new               = sqrt((x_list - circle_x).^2 + (y_list - circle_y).^2) - circle_radius;        
                sum_residuals_new   = sum(R_new.^2);
    
                change = (sum_residuals_new - sum_residuals)/sum_residuals;
                
                % Alternatively, if the relative step size is low it is also said to have converged
                relative_step_size = abs(sum_residuals_new - sum_residuals) / sum_residuals_init;
                            
                if abs(change) < min_step_threshold
                    convergence     = true;
                elseif relative_step_size < min_step_threshold
                    convergence     = true;
                end            
            end
                        
    %% Pratt's least-squares approach for circle fitting %%
        elseif strcmp(Method, 'Pratt')
            % Pratt's method is applied
            [circle_centre, circle_radius] = Pratt_Circle_Fit(point_cloud_matrix, min_step_threshold, max_iterations, weight_list);      
        end

    %% Plot %%
        if Plot == true
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])   
            
            hold on
            grid on

            % Point cloud
            scatter(x_list, y_list, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'DisplayName', 'Point cloud');

            % Circle
            number_coord        = 1e2;
            circle_coord_matrix = Circle_Coordinates(circle_centre, circle_radius, [], number_coord);
            plot(circle_coord_matrix(:, 1), circle_coord_matrix(:, 2), 'LineWidth', 2, 'color', 'b', 'DisplayName', sprintf('%s fitted circle', Method));

            % Axes
            xlabel('x [m]')
            ylabel('y [m]')
            axis equal

            % Legend
            legend('show', 'location', 'eastoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            % Pause message
            fprintf('The circle has been fitted using %s least-squares. The figure will close and script end upon a key-press. \n', Method);
            pause();

            close(1);
        end

end