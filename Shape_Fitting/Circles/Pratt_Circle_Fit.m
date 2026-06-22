% This script applies (weighted) Pratt's least squares to fit a circle to the point cloud matrix
% The step threshold dictates the point at which Newton's method converges, and max iterations the maximum number of iterations that is allowed to take place

function [circle_centre, circle_radius] = Pratt_Circle_Fit(point_cloud_matrix, step_threshold, max_iterations, weights_list)

    %% Point cloud parameters %%
        % Split into dimensions
        [x_list, y_list] = Column_Deal(point_cloud_matrix);
            
        % The weighted centroid of the given coordinates
        x_list_w = x_list .* weights_list;
        y_list_w = y_list .* weights_list;
        
        centroid_x_w = sum(x_list_w) / sum(weights_list);
        centroid_y_w = sum(y_list_w) / sum(weights_list);
    
        % The coordinates are centered
        x_list_c = (x_list_w - centroid_x_w);
        y_list_c = (y_list_w - centroid_y_w);
    
        % The effective radii (z coordinate) at each point is computed
        z_list_c = x_list_c.^2 + y_list_c.^2;
    
        % The moments the points produce around it in 3D    
        M_xx = sum(x_list_c.^2 .* weights_list) / sum(weights_list);
        M_xy = sum(x_list_c .* y_list_c .* weights_list) / sum(weights_list);
        M_xz = sum(x_list_c .* z_list_c .* weights_list) / sum(weights_list);
    
        M_yy = sum(y_list_c.^2 .* weights_list) / sum(weights_list);
        M_yx = sum(y_list_c .* x_list_c .* weights_list) / sum(weights_list);
        M_yz = sum(y_list_c .* z_list_c .* weights_list) / sum(weights_list);
    
        M_zz = sum(z_list_c.^2 .* weights_list) / sum(weights_list);
    
        % These are used for the coefficients of the polynomial describing the circle
        M_z             = M_xx + M_yy;
        Covariance_xy   = M_xx * M_yy - M_xy * M_yx;
    
        A = M_xz^2*M_yy + M_yz^2*M_xx - M_zz*Covariance_xy - 2*M_xz*M_yz*M_xy + M_z*M_z*Covariance_xy;
        B = M_zz*M_z + 4*Covariance_xy*M_z - M_xz^2 - M_yz^2 - M_z^3;
        C = 4*Covariance_xy - 3*M_z^2 - M_zz;

    %% Newton's method %%
        % The loop is initialised
        iter = 0;
    
        D_new = 0;
        E_new = Inf;
    
        while iter <= max_iterations
            iter = iter + 1;
    
            E_old = E_new;
            E_new = A + D_new*(B + D_new*(C + 4*D_new^2));
    
            % A message if the method fails to converge
            if abs(E_new) > abs(E_old)
                error('Convergence cannot be obtained');
            end
    
            F = B + D_new * (C^2 + 16*D_new^2);
            D_old = D_new;
            D_new = D_old - E_new / F;
    
            % Convergence check
            step = (D_new - D_old) / D_old;
    
            if abs(step) < step_threshold
                break
            end
        end
    
        % Check whether the method succeeded
        if iter == max_iterations
            warning('Newton''s method failed to converge');
        end
    
        % The circle parameters
        Determinant     = D_new^2 - D_new * M_z + Covariance_xy;
    
        circle_x_c_w    = (M_xz*(M_yy - D_new) - M_yz*M_xy) / (2*Determinant);
        circle_x        = circle_x_c_w + centroid_x_w;
        circle_y_c_w    = (M_yz*(M_xx - D_new) - M_xz*M_xy) / (2*Determinant);
        circle_y        = circle_y_c_w + centroid_y_w;
    
        circle_centre   = [circle_x, circle_y];
        circle_radius   = sqrt(sum([circle_x_c_w, circle_y_c_w].*[circle_x_c_w, circle_y_c_w]) + M_z + 2*D_new);
end
