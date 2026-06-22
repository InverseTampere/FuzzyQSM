% This script projects the 3D normal distribution onto the plane normal vector (assumed 3rd vector in vector basis) and the plane defined by said vector basis and point
% The projections are given for both the plane and vector, and then the original or plane-based coordinate systems can be chosen (i.e. Projected_Distributions.Vector.Original.mu)

% Note: The coefficients are returned as lists whilst the distribution properties are returned in a cell array
% Note: If the last distribution axis is not parallel to the plane normal vector, it is returned as the last axis of the planar projected distributions

function Projected_Distributions = Multivariate_Normal_Plane_Projection(plane_vector_basis, plane_point, Point_Cloud_Distributions, Plot)

    %% Structure inputs %%
        % Distributions
        distribution_mu_cell        = Point_Cloud_Distributions.distribution_mu_cell;
        distribution_sigmae_cell    = Point_Cloud_Distributions.distribution_sigmae_cell;
        distribution_Sigma_cell     = Point_Cloud_Distributions.distribution_Sigma_cell;
        distribution_axes_cell      = Point_Cloud_Distributions.distribution_axes_cell;
        number_distributions        = Point_Cloud_Distributions.number_distributions;

    %% Manual inputs %%
        % The projection of distributions can be shown in detail one at a time
        Diagnostics             = false ;       % [true, false] 

        % Margin used to check for errors
        diagonality_num_margin  = 1e-6;         % [-] Small value to check if the covariance matrix is diagonal

    %% Projection %%
        %--% The local function is applied to all distributions using cellfun %--%
        Gaussian_Projection_fun  = @(distribution_mu, distribution_axes, distribution_sigmae, distribution_Sigma) Gaussian_Projection(distribution_mu, distribution_axes, distribution_sigmae, distribution_Sigma, plane_vector_basis, plane_point, diagonality_num_margin, Diagnostics);
        
        [plane_coefficients_cell, mu_plane_cell, sigmae_plane_cell, Sigma_plane_cell, distr_axes_plane_cell, mu_plane_3D_cell, Sigma_plane_3D_cell, distr_axes_plane_3D_cell, ...
        vector_coefficients_cell, sigma_vector_cell, mu_vector_cell, mu_vector_3D_cell] = cellfun(Gaussian_Projection_fun, distribution_mu_cell, distribution_axes_cell, distribution_sigmae_cell, distribution_Sigma_cell, 'UniformOutput', false);
    
        %--% The data is collected into one structure %--%
        % Plane distributions
        plane_coefficients_matrix   = vertcat(plane_coefficients_cell{:});
        Plane_Coefficients          = struct('description', '2D ellipse according to alpha*x^2 + beta*y^2 + gamma*xy + delta*x + epsilon*y + omega', 'alpha', plane_coefficients_matrix(:, 1), 'beta', plane_coefficients_matrix(:, 2), 'gamma', plane_coefficients_matrix(:, 3), 'delta', plane_coefficients_matrix(:, 4), 'epsilon', plane_coefficients_matrix(:, 5), 'omega', plane_coefficients_matrix(:, 6));

        Proj_Plane_Distributions    = struct('mu', {mu_plane_cell}, 'Sigma', {Sigma_plane_cell}, 'sigmae', {sigmae_plane_cell}, 'distr_axes', {distr_axes_plane_cell});
        Orig_Plane_Distributions    = struct('mu', {mu_plane_3D_cell}, 'Sigma', {Sigma_plane_3D_cell}, 'sigmae', {sigmae_plane_cell}, 'distr_axes', {distr_axes_plane_3D_cell});

        Plane_Distributions         = struct('Coefficients', Plane_Coefficients, 'Projection', Proj_Plane_Distributions, 'Original', Orig_Plane_Distributions);

        % Vector distributions
        vector_coefficients_matrix  = vertcat(vector_coefficients_cell{:}); 
        Vector_Coefficients         = struct('description', '2D ellipse according to alpha*x^2 + beta*y^2 + gamma*xy + delta*x + epsilon*y + omega, 1D ellipse according to zeta*z^2 + eta*z + theta', 'alpha', vector_coefficients_matrix(:, 1), 'beta', vector_coefficients_matrix(:, 2), 'gamma', vector_coefficients_matrix(:, 3), 'delta', vector_coefficients_matrix(:, 4), 'epsilon', vector_coefficients_matrix(:, 5), 'omega', vector_coefficients_matrix(:, 6), 'zeta', vector_coefficients_matrix(:, 7), 'eta', vector_coefficients_matrix(:, 8), 'theta', vector_coefficients_matrix(:, 9));
        
        Proj_Vector_Distributions   = struct('mu', {mu_vector_cell}, 'sigma', {sigma_vector_cell});
        Orig_Vector_Distributions   = struct('mu', {mu_vector_3D_cell}, 'sigma', {sigma_vector_cell});

        Vector_Distributions        = struct('Coefficients', Vector_Coefficients, 'Projection', Proj_Vector_Distributions, 'Original', Orig_Vector_Distributions);

        % The integrated distribution properties are put in one structure for the vector and plane
        Projected_Distributions     = struct('vector_basis', plane_vector_basis, 'plane_point', plane_point, 'number_distributions', number_distributions, 'Plane', Plane_Distributions, 'Vector', Vector_Distributions);
        
    %% Plot %%
        if Plot == true
            num_coord       = 1e2;          % [-] Number of coordinates for each distribution
            num_STD         = 1.00;         % [-] A smaller factor can make it easier to distinguish distributions, but more difficult to see their shapes  

            % Used colours
            orig_3D_colour      = 'k';      % Original 3D distributions
            plane_proj_colour   = 'b';      % Projections onto the plane
            vector_proj_colour  = 'r';      % Projections onto the vector

            %--% Original coordinate frame %--%
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            grid on
            hold on
            
            % 3D Distributions
            for d = 1: number_distributions
                % 3D distribution properties
                distr_mu        = distribution_mu_cell{d};
                distr_radii     = num_STD * distribution_sigmae_cell{d};
                distr_axes      = distribution_axes_cell{d};
                
                % Ellipsoidal surface
                [ellipsoid_coord_matrix, num_coord] = Ellipsoid_Coordinate_Generator(distr_mu, distr_radii, distr_axes, num_coord);

                x_ellipsoid = reshape(ellipsoid_coord_matrix(:, 1), sqrt(num_coord) * [1, 1]);
                y_ellipsoid = reshape(ellipsoid_coord_matrix(:, 2), sqrt(num_coord) * [1, 1]);
                z_ellipsoid = reshape(ellipsoid_coord_matrix(:, 3), sqrt(num_coord) * [1, 1]);

                surf_distr = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', orig_3D_colour, 'FaceAlpha', 0.25, 'DisplayName', sprintf('Distribution, %.3g %s', num_STD, '\sigma'));

                % mu
                sc_mu = scatter3(distr_mu(1), distr_mu(2), distr_mu(3), 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', orig_3D_colour, 'DisplayName', '\mu');

                if d > 1
                    surf_distr.HandleVisibility = 'Off';
                    sc_mu.HandleVisibility      = 'Off';
                end
            end
            
            % Projection vector basis
            mu_plane_3D_matrix          = vertcat(mu_plane_3D_cell{:});
            delta_mu_plane_3D_matrix    = mu_plane_3D_matrix - plane_point;
            delta_mu_plane_3D_norm_list = sqrt(sum(delta_mu_plane_3D_matrix.^2, 2));
            plane_extent                = mean(delta_mu_plane_3D_norm_list);

            mu_vector_3D_matrix             = vertcat(mu_vector_3D_cell{:});
            delta_mu_vector_3D_matrix       = mu_vector_3D_matrix - plane_point;
            delta_mu_vector_3D_norm_list    = sqrt(sum(delta_mu_vector_3D_matrix.^2, 2));
            orth_vector_extent              = mean(delta_mu_vector_3D_norm_list);

            [a_vec, b_vec, c_vec]   = deal(plane_vector_basis(1, :), plane_vector_basis(2, :), plane_vector_basis(3, :));
            
            plot3(plane_point(1) + plane_extent * [0, a_vec(1)], plane_point(2) + plane_extent * [0, a_vec(2)], plane_point(3) + plane_extent * [0, a_vec(3)], 'LineWidth', 2, 'LineStyle', ':', 'color', plane_proj_colour, 'DisplayName', 'Plane vector a');
            plot3(plane_point(1) + plane_extent * [0, b_vec(1)], plane_point(2) + plane_extent * [0, b_vec(2)], plane_point(3) + plane_extent * [0, b_vec(3)], 'LineWidth', 2, 'LineStyle', ':', 'color', plane_proj_colour, 'DisplayName', 'Plane vector b');
            plot3(plane_point(1) + orth_vector_extent * [0, c_vec(1)], plane_point(2) + orth_vector_extent * [0, c_vec(2)], plane_point(3) + orth_vector_extent * [0, c_vec(3)], 'LineWidth', 2, 'LineStyle', ':', 'color', vector_proj_colour, 'DisplayName', 'Orth. vector c');

            % Plane projected distributions
            for d = 1 : number_distributions
                % 3D plane projected distribution properties
                mu_plane            = mu_plane_3D_cell{d};
                ellipse_radii_plane = num_STD * sigmae_plane_cell{d};
                distr_axes_plane    = distr_axes_plane_3D_cell{d};

                % Elliptical distribution
                ellipse_coord_matrix    = Ellipse_Coordinate_Generator(mu_plane, distr_axes_plane, ellipse_radii_plane, num_coord);
                pl_distr_plane          = plot3(ellipse_coord_matrix(:, 1), ellipse_coord_matrix(:, 2), ellipse_coord_matrix(:, 3), 'LineWidth', 2, 'color', plane_proj_colour, 'DisplayName', sprintf('Plane proj. distr., %.3g %s', num_STD, '\sigma'));

                % mu
                sc_mu_plane = scatter3(mu_plane(1), mu_plane(2), mu_plane(3), 'filled', 'MarkerFaceColor', plane_proj_colour, 'MarkerEdgeColor', 'none', 'DisplayName', 'Plane proj. \mu');

                if d > 1
                    pl_distr_plane.HandleVisibility = 'Off';
                    sc_mu_plane.HandleVisibility    = 'Off';
                end
            end

            % Vector projected distributions
            for d = 1 : number_distributions
                % 1D vector projected distribution properties
                mu_vector       = mu_vector_3D_cell{d};
                extent_vector   = sigma_vector_cell{d};

                % Extent on the vector
                pl_distr_vec = plot3(mu_vector(1) + extent_vector*c_vec(1)*[-1, 1], mu_vector(2) + extent_vector*c_vec(2)*[-1, 1], mu_vector(3) + extent_vector*c_vec(3)*[-1, 1], 'LineWidth', 2, 'color', vector_proj_colour, 'DisplayName', sprintf('Vector proj. distr., %.3g %s', num_STD, '\sigma'));

                % mu
                sc_mu_vec = scatter3(mu_vector(1), mu_vector(2), mu_vector(3), 'filled', 'MarkerFaceColor', vector_proj_colour, 'MarkerEdgeColor', 'none', 'DisplayName', 'Vector proj. \mu');

                if d > 1
                    pl_distr_vec.HandleVisibility   = 'Off';
                    sc_mu_vec.HandleVisibility      = 'Off';
                end
            end

            % Axes
            xlabel('x [m]');
            ylabel('y [m]');
            zlabel('z [m]');
            axis equal

            view(45, 45);

            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    
            
            %--% Plane projected coordinate frame %--%
            figure(2)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            grid on
            hold on

            for d = 1 : number_distributions
                % 2D plane projected distribution properties
                mu_plane            = mu_plane_cell{d};
                ellipse_radii_plane = num_STD * sigmae_plane_cell{d};
                distr_axes_plane    = distr_axes_plane_cell{d};

                % Elliptical distribution
                ellipse_coord_matrix    = Ellipse_Coordinate_Generator(mu_plane, distr_axes_plane, ellipse_radii_plane, num_coord);
                pl_distr_plane          = plot(ellipse_coord_matrix(:, 1), ellipse_coord_matrix(:, 2), 'LineWidth', 2, 'color', plane_proj_colour, 'DisplayName', sprintf('Plane proj. distr., %.3g %s', num_STD, '\sigma'));

                % mu
                sc_mu_plane = scatter(mu_plane(1), mu_plane(2), 'filled', 'MarkerFaceColor', plane_proj_colour, 'MarkerEdgeColor', 'none', 'DisplayName', 'Plane proj. \mu');

                if d > 1
                    pl_distr_plane.HandleVisibility = 'Off';
                    sc_mu_plane.HandleVisibility    = 'Off';
                end            
            end

            % Axes
            xlabel('a [m]');
            ylabel('b [m]');

            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off  

            % Pause message
            disp('Projection of the MVNs has finished. The script will end and figures will close upon a key-press.');
            pause();
            close(1:2);
        end

    %% Local projection function %%
    function  [plane_coefficients, mu_plane, sigmae_plane, Sigma_plane, distr_axes_plane, mu_plane_3D, Sigma_plane_3D, distr_axes_plane_3D, ...
               vector_coefficients, sigma_vector, mu_vector, mu_vector_3D] = Gaussian_Projection(distribution_mu, distribution_axes, distribution_sigmae, distribution_Sigma, plane_vector_basis, plane_point, diagonality_num_margin, Diagnostics)
           
        %% Transformation to the plane's coordinate frame %%
            % Expected value vector
            distribution_mu_t   = distribution_mu - plane_point;
            mu                  = (plane_vector_basis * distribution_mu_t')';

            % Covariance matrix
            Sigma               = plane_vector_basis * distribution_Sigma * plane_vector_basis';

        %% 3D ellipsoid parametrisation %%
            % The exponential of a Gaussian is -1/2 times the Mahalanobis distance (M) squared, which can be rewritten to easily parametrise a 3D ellipsoid where
            % -1/2 M^2 = -1/2 (Pxx*x^2 + Pyy*y^2 + Pzz*z^2) + Ax + By + Cz - Pxy*xy - Pxz*xy - Pyz*yz + D

            % P is the precision matrix, i.e. the inverse of the covariance matrix
            P       = Inverse_of_3x3_Matrix(Sigma);
            num_dim = size(P, 1);

            % Together with the expected value it leads to the constants A through D
            [i_vec, j_vec, k_vec]   = deal(P(1, :), P(2, :), P(3, :));
            [A, B, C]               = deal(dot(mu, i_vec), dot(mu, j_vec), dot(mu, k_vec));
            D                       = -1/2 * dot(diag(P), mu.^2) - P(1, 2)*mu(1)*mu(2) - P(1, 3)*mu(1)*mu(3) - P(2, 3)*mu(2)*mu(3);

        %% Integration onto plane %%
            % Integration onto the plane is a single integral in the third dimension for the transformed distribution using the following set of coefficients    
            % The coefficients then dictate the 2D ellipse according to alpha*x^2 + beta*y^2 + gamma*xy + delta*x + epsilon*y + omega
            alpha_z     = 1/2*(P(1, 3)^2/P(3, 3) - P(1, 1));
            beta_z      = 1/2*(P(2, 3)^2/P(3, 3) - P(2, 2));
            gamma_z     = P(1, 3)*P(2, 3)/P(3, 3) - P(1, 2);            % Measure of the diagonality. Zero if P is purely diagonal
            delta_z     = A - C*P(1, 3)/P(3, 3);
            epsilon_z   = B - C*P(2, 3)/P(3, 3);
            omega_z     = C^2/(2*P(3, 3)) + D;     
            
            plane_coefficients = [alpha_z, beta_z, gamma_z, delta_z, epsilon_z, omega_z];
           
            if abs(gamma_z) < diagonality_num_margin                    % Margin on the diagonality in case of numerical rounding             
                % If the covariance matrix is diagonal, planar projection is simply the first two dimensions        
                mu_plane            = mu(1 : num_dim - 1);
                Sigma_plane         = Sigma(1 : num_dim - 1, 1 : num_dim - 1);
                sigmae_plane        = sqrt(diag(Sigma_plane))';
                distr_axes_plane    = eye(num_dim - 1);
            else
                % Integrated distribution properties on the plane
                mu_plane_x      = (epsilon_z*gamma_z - 2*beta_z*delta_z) / (4*alpha_z*beta_z - gamma_z^2);
                mu_plane_y      = -(delta_z + 2*alpha_z*mu_plane_x) / gamma_z;
                mu_plane        = [mu_plane_x, mu_plane_y];
    
                psi             = 1/(gamma_z^2 - 4*alpha_z*beta_z);
                Sigma_plane     = psi * [2*beta_z, -gamma_z; -gamma_z, 2*alpha_z];
    
                % The axes are the eigenvectors of the covariance matrix, and the standard deviations the square roots of the eigenvalues
                lambda_list         = psi*(alpha_z + beta_z + [-1; 1]*sqrt(alpha_z^2 + beta_z^2 + gamma_z^2 - 2*alpha_z*beta_z));        
                sigmae_plane        = sqrt(lambda_list)';
    
                distr_axes_plane    = 1./sqrt(psi^2*gamma_z^2 + (2*psi*beta_z - lambda_list).^2) .* [repmat(psi * gamma_z, [2, 1]), 2*psi*beta_z - lambda_list];

                % The axis with the largest x-component is changed to be first
                [~, first_axis_ind] = max(abs(distr_axes_plane(:, 1)));
                second_axis_ind     = setdiff(1:num_dim - 1, first_axis_ind);

                distr_axes_plane    = distr_axes_plane([first_axis_ind, second_axis_ind], :);
                sigmae_plane        = sigmae_plane([first_axis_ind, second_axis_ind]);
            end

            % Transformed to the original coordinate frame        
            mu_plane_3D_r       = [mu_plane, 0];
            mu_plane_3D_t       = (plane_vector_basis' * mu_plane_3D_r')';
            mu_plane_3D         = mu_plane_3D_t + plane_point;

            Sigma_plane_3D_r    = [[Sigma_plane, [0; 0]]; [0, 0, 0]];
            Sigma_plane_3D      = plane_vector_basis' * Sigma_plane_3D_r * plane_vector_basis;

            distr_axes_plane_3D = (plane_vector_basis' * [[distr_axes_plane; zeros(1, 2)], [0; 0; 1]]')';
            distr_axes_plane_3D = distr_axes_plane_3D(1 : num_dim - 1, :);                                  % The last is the plane normal vector and is removed

        %% Integration onto vector %%
            % Integration onto the vector is a double integral in the first and second dimension for the transformed distribution using the following set of coefficients    
            % The first integral is performed with respect to x and produces an ellipse on the y-z plane with the following coefficients
            alpha_x     = 1/2*(P(1, 2)^2/P(1, 1) - P(2, 2));
            beta_x      = 1/2*(P(1, 3)^2/P(1, 1) - P(3, 3));
            gamma_x     = P(1, 2)*P(1, 3)/P(1, 1) - P(2, 3);
            delta_x     = B - A*P(1, 2)/P(1, 1);
            epsilon_x   = C - A*P(1, 3)/P(1, 1);
            omega_x     = A^2/(2*P(1, 1)) + D;    

            % The second integral uses coefficients which are a function of the first set forming a '1D ellipse' zeta*z^2 + eta*z + theta
            zeta_y      = beta_x - gamma_x^2 / (4*alpha_x);
            eta_y       = epsilon_x - delta_x*gamma_x / (2*alpha_x);
            theta_y     = omega_x - delta_x^2 / (4*alpha_x);
            
            vector_coefficients = [alpha_x, beta_x, gamma_x, delta_x, epsilon_x, omega_x, zeta_y, eta_y, theta_y];

            % The twice-integrated distribution's properties
            sigma_vector    = sqrt(-1/(2*zeta_y));
            mu_vector       = -eta_y / (2*zeta_y);

            % The expected value is transformed to the original coordinate frame
            mu_vector_3D_r  = [0, 0, mu_vector];
            mu_vector_3D_t  = (plane_vector_basis' * mu_vector_3D_r')';
            mu_vector_3D    = mu_vector_3D_t + plane_point;

        %% Diagnostics plot %%
            if Diagnostics == true
                % Discretisation
                number_samples  = 1e3;        
                m_STD           = 1;

                %--% Original coordinate frame %--%
                figure(1)
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])    

                hold on
                grid on

                % The distribution
                distribution_radii = distribution_sigmae * m_STD;
                [distr_coord_matrix, number_coord] = Ellipsoid_Coordinate_Generator(distribution_mu, distribution_radii, distribution_axes, number_samples);

                a_distribution  = reshape(distr_coord_matrix(:, 1), sqrt(number_coord) * [1, 1]);
                b_distribution  = reshape(distr_coord_matrix(:, 2), sqrt(number_coord) * [1, 1]);
                c_distribution  = reshape(distr_coord_matrix(:, 3), sqrt(number_coord) * [1, 1]);

                surf_distr_name = sprintf('Distribution, %.2g %s', m_STD, '\sigma');
                surf(a_distribution, b_distribution, c_distribution, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'DisplayName', surf_distr_name);

                % Projected plane distributions
                radii_plane         = m_STD * sigmae_plane;
                plane_distr_coord   = Ellipse_Coordinate_Generator(mu_plane_3D, distr_axes_plane_3D, radii_plane, number_samples);

                pl_plane_name = sprintf('Plane distribution, %.2g %s', m_STD, '\sigma');
                plot3(plane_distr_coord(:, 1), plane_distr_coord(:, 2), plane_distr_coord(:, 3), 'LineWidth', 2, 'color', 'g', 'DisplayName', pl_plane_name);

                % The plane
                plane_normal_vector     = plane_vector_basis(num_dim, :);
                plane_corner_matrix     = Plane_Corner_Points(plane_normal_vector, plane_point, [distr_coord_matrix; plane_distr_coord]);
                pa_plane_name           = 'Plane';
                patch(plane_corner_matrix(:, 1), plane_corner_matrix(:, 2), plane_corner_matrix(:, 3), 'g', 'FaceAlpha', 0.25, 'DisplayName', pa_plane_name);   

                % The plane normal vector
                vector_length   = 2*max(distribution_radii);        
                pl_vector_name  = 'Normal vector';
                plot3(plane_point(1) + vector_length * [-1, 1] * plane_normal_vector(1), plane_point(2) + vector_length * [-1, 1] * plane_normal_vector(2), plane_point(3) + vector_length * [-1, 1] * plane_normal_vector(3), 'LineWidth', 2, 'color', 'm', 'DisplayName', pl_vector_name);

                % Vector distribution expected value
                sc_mu_vec_name = '\mu_{v}';
                scatter3(mu_vector_3D(1), mu_vector_3D(2), mu_vector_3D(3), 'filled', 'MarkerFaceColor', 'm', 'DisplayName', sc_mu_vec_name);

                % Axes
                xlabel('a [m]');
                ylabel('b [m]');
                zlabel('c [m]');

                set(gca, 'DataAspectRatio', [1, 1, 1]);

                % Text
                legend('show', 'location', 'northoutside');

                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                % The perspective
                view(45, 45);

                %--% Plane-aligned coordinate frame %--%
                figure(2)
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])    

                hold on
                grid on

                % The distribution
                distribution_radii      = distribution_sigmae * m_STD;
                distribution_axes_xyz   = (plane_vector_basis * distribution_axes')';
                [distr_coord_matrix, number_coord] = Ellipsoid_Coordinate_Generator(mu, distribution_radii, distribution_axes_xyz, number_samples);

                x_distribution  = reshape(distr_coord_matrix(:, 1), sqrt(number_coord) * [1, 1]);
                y_distribution  = reshape(distr_coord_matrix(:, 2), sqrt(number_coord) * [1, 1]);
                z_distribution  = reshape(distr_coord_matrix(:, 3), sqrt(number_coord) * [1, 1]);

                surf_distr_name = sprintf('Distribution, %.2g %s', m_STD, '\sigma');
                surf(x_distribution, y_distribution, z_distribution, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'DisplayName', surf_distr_name);

                % Last vector
                last_vector         = distribution_axes_xyz(num_dim, :);
                last_vector_scaled  = distribution_radii(num_dim) * last_vector; 

                plot3(mu(1) + [0, last_vector_scaled(1)], mu(2) + [0, last_vector_scaled(2)], mu(3) + [0, last_vector_scaled(3)], 'LineWidth', 2, 'color', 'b', 'LineStyle', ':', 'DisplayName', 'Last vector');

                % Projected plane distributions
                radii_plane         = m_STD * sigmae_plane;
                plane_distr_coord   = Ellipse_Coordinate_Generator(mu_plane, distr_axes_plane, radii_plane, number_samples);

                pl_plane_name   = sprintf('Plane distribution, %.2g %s', m_STD, '\sigma');
                plot3(plane_distr_coord(:, 1), plane_distr_coord(:, 2), zeros(number_samples, 1), 'LineWidth', 2, 'color', 'g', 'DisplayName', pl_plane_name);

                % Vector distribution expected value
                sc_mu_vec_name      = '\mu_{v}';
                scatter3(0, 0, mu_vector, 'filled', 'MarkerFaceColor', 'm', 'DisplayName', sc_mu_vec_name);

                % Axes
                xlabel('x [m]');
                ylabel('y [m]');
                zlabel('z [m]');

                set(gca, 'DataAspectRatio', [1, 1, 1]);

                % Text
                legend('show', 'location', 'northoutside');

                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                % The perspective
                view(45, 45);

                %--% Pause %--%
                disp('The script will continue and the figures will close upon a button-press');
                pause();

                close(1:2); 
            end
    end
end


