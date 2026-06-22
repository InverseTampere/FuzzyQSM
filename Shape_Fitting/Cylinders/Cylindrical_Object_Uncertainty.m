% This script computes the radial and range uncertainty of laser beams hitting a cylindrical object

% The incidence angle (and thus range uncertainty) can grow to infinity for grazing hits
% As such, for points for which the maximum incidence angle is exceeded the uncertainty is substituted
% max_incidence_angle can bet set to >pi/2 if this is undesired

function [sigma_radial_cell, sigma_prop_cell, incidence_angle_cell, beam_range_cell, beamwidth_cell] = Cylindrical_Object_Uncertainty(cyl_centre, cyl_direction, Scanner_loc_cell, Scanner_Parameters, Point_Cloud_Coord)

    %% Inputs %%
        % Scanner parameters
        max_incidence_angle         = Scanner_Parameters.max_incidence_angle;     

        % The point cloud
        point_cloud_cell            = Point_Cloud_Coord.point_cloud_cell;

    %% The point cloud uncertainty %%
        % The incidence angles
        number_scanners         = length(Scanner_loc_cell);
        incidence_angle_cell    = cell(1, number_scanners);
        
        for s = 1 : number_scanners
            % The data for this scanner 
            point_matrix        = point_cloud_cell{s};
            number_points       = size(point_matrix, 1);
            
            % The beam vectors
            Scanner_loc         = Scanner_loc_cell{s};
            scanner_loc_matrix  = repmat(Scanner_loc, number_points, 1);
            beam_vector_matrix  = point_matrix - scanner_loc_matrix;
            
            beam_range_list     = sqrt(sum(beam_vector_matrix.^2, 2));
            beam_vector_matrix  = beam_vector_matrix ./ beam_range_list;
                        
            %--% Incidence angles %--%
            % Projection of the points onto the cylinder axis
            [proj_point_matrix, ~, omega_list]  = Point_to_Vector_Projection(point_matrix, cyl_direction, cyl_centre);
            proj_point_vector_matrix            = (point_matrix - proj_point_matrix) ./ omega_list;
            
            % The angles between the laser beam and the vector from the point to the cylinder axis
            beta_angle_list         = acos(dot(proj_point_vector_matrix, beam_vector_matrix, 2));
            incidence_angle_list    = min(beta_angle_list, pi - beta_angle_list);
                
            % To ensure that excessively large uncertainties are not present, the incidence angle is capped
            incidence_angle_list    = min(incidence_angle_list, max_incidence_angle);          
            incidence_angle_cell{s} = incidence_angle_list;
        end

        % The uncertainties
        [sigma_radial_cell, sigma_prop_cell, beam_range_cell, beamwidth_cell] = Gaussian_Beam_Uncertainty(Point_Cloud_Coord, Scanner_Parameters, Scanner_loc_cell, incidence_angle_cell);

end