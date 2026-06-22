% The distributions are projected onto the cross-section with the (approximate) propagation axis being -y

function [proj_mu_matrix, proj_sigmae_matrix] = Distribution_Cylinder_Cross_Section_Projection(Point_Cloud_Distributions, Scanner_loc_cell, cylinder_centre, cylinder_direction)

    %% Structure inputs %%
        number_distributions_list   = Point_Cloud_Distributions.number_distributions_list;

    %% Manual inputs %%
        Proj_Diagnostics            = false;            % [true, false] Projection onto the cross-sectional coordinate frame

    %% Projection %%
        % As the propagation axis is diferent for each scanner, they each have their own vector basis
        proj_vector_basis_cell = Vector_Basis_Cylinder_Cross_Section_Projection(cylinder_centre, cylinder_direction, Scanner_loc_cell);
            
        % The transformed properties are stored in matrix format
        cumul_beams_list        = [0, cumsum(number_distributions_list)];
        total_num_distributions = sum(number_distributions_list);
        num_dim                 = length(cylinder_centre);
        number_scanners         = length(Scanner_loc_cell);

        proj_mu_matrix      = zeros(total_num_distributions, num_dim - 1);
        proj_sigmae_matrix  = zeros(total_num_distributions, num_dim - 1);

        for s = 1 : number_scanners
            % This scanner's indices
            ind_first   = cumul_beams_list(s) + 1;
            ind_last    = cumul_beams_list(s + 1);

            % The vector basis
            projection_vector_basis = proj_vector_basis_cell{s};       

            % Distributions projected on the cross-sectional plane
            Projected_Distributions_scanner = Multivariate_Normal_Plane_Projection(projection_vector_basis, cylinder_centre, Point_Cloud_Distributions, Proj_Diagnostics);

            % Projected mu and sigmae
            mu_Q_cell                               = Projected_Distributions_scanner.Plane.Projection.mu;
            mu_Q_matrix                             = vertcat(mu_Q_cell{:});
            proj_mu_matrix(ind_first : ind_last, :) = mu_Q_matrix(ind_first : ind_last, :);

            sigmae_Q_cell                               = Projected_Distributions_scanner.Plane.Projection.sigmae;
            sigmae_Q_matrix                             = vertcat(sigmae_Q_cell{:});
            proj_sigmae_matrix(ind_first : ind_last, :) = sigmae_Q_matrix(ind_first : ind_last, :);
        end
end