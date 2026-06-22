% This script generates the corner points of a plane in 3D
% The extent is dictated by the data matrix

function plane_corner_matrix = Plane_Corner_Points(plane_normal_vector, plane_point, data_matrix)
        
    %% Manual inputs %%
        bound_margin = 1.2;     % Factor by which the bounds are extended

    %% Bounds of the data %%
        UB_list     = max([data_matrix; plane_point]);           % Note that the plane point is included  
        LB_list     = min([data_matrix; plane_point]);
        
        ampl_list   = UB_list - LB_list;
        
        UB_list     = UB_list  + (bound_margin - 1)/2 * ampl_list;
        LB_list     = LB_list  - (bound_margin - 1)/2 * ampl_list;

    %% Determine the 'locked' dimension %%
        % The first nonzero component of the plane normal vector
        plane_normal_vector = plane_normal_vector / norm(plane_normal_vector);
        nonzero_component   = find(abs(plane_normal_vector) == max(abs(plane_normal_vector)), 1);

    %% Corner coordinates (x-y plane) %%
        num_dim     = length(plane_normal_vector);
        num_corners = 4;
        
        plane_corner_matrix = zeros(num_corners, num_dim);
        
        % The free dimensions are given the bounds as corner points
        free_dimensions = setdiff(1 : num_dim, nonzero_component);
        
        for d = 1 : num_dim - 1
            dimension = free_dimensions(d);
            
            if d == 1
                corner_list = [UB_list(dimension), LB_list(dimension), LB_list(dimension), UB_list(dimension)];
            else
                corner_list = [UB_list(dimension), UB_list(dimension), LB_list(dimension), LB_list(dimension)];
            end
            
            plane_corner_matrix(:, dimension) = corner_list;
        end
    
        % The locked dimension's coordinates depend on those of the free dimensions
        plane_normal_vector = plane_normal_vector / norm(plane_normal_vector);
        plane_constant      = dot(plane_normal_vector, plane_point);
        plane_normal_matrix = repmat(plane_normal_vector, [num_corners, 1]);     % The normal vector is repeated for the dot product
        
        locked_list = (plane_constant - dot(plane_normal_matrix, plane_corner_matrix, 2)) / plane_normal_vector(nonzero_component);     % Note that the locked dimension is effectively disregarded in the dot product, as it is 0 so far
        
        plane_corner_matrix(:, nonzero_component) = locked_list;
        
end