% This script computes the inverse of a 3x3 matrix analytically

function inverse_matrix = Inverse_of_3x3_Matrix(matrix)

    % Each entry in the matrix is assigned a variable
    [a, b, c] = Column_Deal(matrix(1, :));
    [d, e, f] = Column_Deal(matrix(2, :));
    [g, h, i] = Column_Deal(matrix(3, :));

    % The determinant and adjugate
    determinant = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
    adjugate    = [e*i - f*h, c*h - b*i, b*f - c*e; ...
                   f*g - d*i, a*i - c*g, c*d - a*f; ...
                   d*h - e*g, b*g - a*h, a*e - b*d];

    % The inverse matrix
    inverse_matrix = adjugate / determinant;
end