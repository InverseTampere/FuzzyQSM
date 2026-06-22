% Two vectors are given which are the first, second or third vector of a 3D vector basis
% The remaining vector is determined using the cross product

function third_vector = Cross_Product_3D(first_vector, first_dim, second_vector, second_dim)

    % Order of cross products
    cross_product_matrix = repmat(1 : 3, [3, 1]) + (0:2)';      % The dimension goes up by one every time
    cross_product_matrix = mod(cross_product_matrix, 3);        % There are only three dimensions
    cross_product_matrix(cross_product_matrix == 0) = 3;        % The modulus of 3 with 3 is zero, but replaced with 3 here

    % The remaining vector dimension
    third_dim = setdiff(1 : 3, [first_dim, second_dim]);

    % Associated row in the cross product matrix
    product_bool    = cross_product_matrix(:, 3) == third_dim;
    cross_product   = cross_product_matrix(product_bool, :);

    % Performing the cross product
    first_vector    = first_vector / norm(first_vector);        % Ensured to be of unit length
    second_vector   = second_vector / norm(second_vector);      

    if first_dim == cross_product(1)
        third_vector = cross(first_vector, second_vector);
    else
        third_vector = cross(second_vector, first_vector);
    end
end