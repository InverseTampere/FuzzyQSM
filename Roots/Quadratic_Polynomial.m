% The two roots to the quadratic polynomial ax^2 + bx + c are found using the simple equation
% The complex output gives all roots regardless of realness

function [x_roots_real, number_real_roots, x_roots_complex] = Quadratic_Polynomial(a, b, c)
    
    %% Roots %%
        % if a is non-zero, the quadratic polynomial is solved
        if a ~= 0
            % Discriminant
            D = b^2 - 4*a*c;
            
            % The realness and uniqueness of the roots follows from the sign of the discriminant
            if D == 0       % If it is zero, there is 1 real root
                x_roots_real        = -b/(2*a);
                x_roots_complex     = x_roots_real;
                
            elseif D > 0    % If it is positive, there are 2 real roots
                x_roots_real        = (-b + [-1; 1] * sqrt(b^2 - 4*a*c)) / (2*a);
                x_roots_complex     = x_roots_real;
                
            else            % If it is negative, there are 2 complex roots
                x_roots_real        = [];
                x_roots_complex     = (-b + [-1; 1] * sqrt(b^2 - 4*a*c)) / (2*a);
            end
    
            number_real_roots = length(x_roots_real);
        
        % As long as b is non-zero it is the linear equation bx + c = 0
        elseif b ~= 0
            number_real_roots   = 1;
            x_roots_real        = -c/b;
            x_roots_complex     = -c/b;

        % Otherwise a solution cannot be found
        else
            error('The quadratic equation cannot be solved as a = b = 0.');
        end
end