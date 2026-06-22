% The beamwidth and resulting Gaussian radial uncertainty are determined here. If the incidence angle is given, the propagation uncertainty is determined as well

function [sigma_radial_cell, sigma_prop_cell, beam_range_cell, beamwidth_cell] = Gaussian_Beam_Uncertainty(Point_Cloud_Coord, Scanner_Parameters, Scanner_loc_cell, incidence_angle_cell)
    
    %% Inputs %%
        % Scanner parameters
        beam_divergence             = Scanner_Parameters.beam_divergence;        
        beam_exit_diameter          = Scanner_Parameters.beam_exit_diameter;          
        sigma_range_device          = Scanner_Parameters.sigma_range_device;       
        
        % The point cloud
        point_cloud_cell            = Point_Cloud_Coord.point_cloud_cell;

    %% Beam properties %%
        % Computed for each scanner
        number_scanners     = length(Scanner_loc_cell);
        sigma_radial_cell   = cell(1, number_scanners);
        sigma_prop_cell     = cell(1, number_scanners);
        beam_range_cell     = cell(1, number_scanners);
        beamwidth_cell      = cell(1, number_scanners);

        for s = 1 : number_scanners
            % This scanner's data
            point_cloud_matrix  = point_cloud_cell{s};
            scanner_loc         = Scanner_loc_cell{s};

            % The beam vectors
            beam_vector_matrix  = point_cloud_matrix - scanner_loc;
            
            beam_range_list     = sqrt(sum(beam_vector_matrix.^2, 2));
            beam_range_cell{s}  = beam_range_list;
            
            beamwidth_list      = beam_exit_diameter + 2*beam_range_list * tan(beam_divergence);
            beamwidth_cell{s}   = beamwidth_list;
            
            % Radial uncertainty
            sigma_radial_list       = beamwidth_list / 4;           % As the beamwidth is defined to be 4 standard deviations, i.e. 2 in each direction
            sigma_radial_cell{s}    = sigma_radial_list;

            % Propagation uncertainty
            if ~isempty(incidence_angle_cell)
                incidence_angle_list    = incidence_angle_cell{s};
                sigma_prop_list         = sigma_range_device + sigma_radial_list .* tan(incidence_angle_list);
            else
                sigma_prop_list         = [];
            end

            sigma_prop_cell{s}      = sigma_prop_list;
        end
end