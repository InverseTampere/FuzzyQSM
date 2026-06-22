% This script computes the coordinates of an ellipse, given its radii, centre and axes

function ellipse_coord_matrix = Ellipse_Coordinate_Generator(ellipse_centre, ellipse_axes, ellipse_radii, number_coord)

    %% The ellipse coordinates %%
        % Discretised angle around the centre
        phi_list = linspace(0, 2*pi, number_coord)';
        
        % Ellipse coordinates
        ellipse_coord_matrix = ellipse_centre + ellipse_radii(1)*ellipse_axes(1, :).*cos(phi_list) + ellipse_radii(2)*ellipse_axes(2, :).*sin(phi_list);
                
end