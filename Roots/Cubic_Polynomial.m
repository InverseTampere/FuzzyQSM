% The roots of the cubic polynomial a*x^3 + b*x^2 + c*x + d are found
% They are split into a real and complex list, where the number of complex roots is always three

function [x_roots_real, number_real_roots, x_roots_complex] = Cubic_Polynomial(a, b, c, d, Plot, Print)

    %% Roots %%
        % If a is nonzero, the cubic polynomial is solved
        if a ~= 0
            % Split equation coefficients
            delta_0 = b^2 - 3*a*c;
            delta_1 = 2*b^3 - 9*a*b*c + 27*a^2*d;
    
            % Sign used in cubic root coefficient
            if delta_0 == 0
                s = -1;
            else
                s = 1;
            end
    
            C = (1/2 * (delta_1 + s*sqrt(delta_1^2 - 4*delta_0^3)))^(1/3);
    
            % Three roots
            k_list  = 0:2;
            xi_list = ((-1 + sqrt(-3)) / 2).^(k_list)';
    
            x_roots_complex = -1/(3*a) * (b + xi_list * C + delta_0 ./ (xi_list * C));
    
            % Complex roots are removed
            tiny_imag_components = abs(imag(x_roots_complex)) < 1e-6;                               % A very minor imaginary component may be present due to Matlab rounding that is removed
            x_roots_complex(tiny_imag_components) = real(x_roots_complex(tiny_imag_components));
            
            x_roots_real        = x_roots_complex(tiny_imag_components);         
            number_real_roots   = length(x_roots_real);

        % Otherwise it is solved as a quadratic polynomial
        else
            [x_roots_real, number_real_roots, x_roots_complex]  = Quadratic_Polynomial(b, c, d);
        end

    %% Printed messages %%
        if Print == true
            fprintf('Roots: (%g are real) \n', number_real_roots);
            
            for r = 1 : 3
                root = x_roots_complex(r);
                fprintf('   %.3g + %.3gi \n', real(root), imag(root));
            end
        end
        
    %% Plot %%
        if Plot == true
            % Amplitude of the variable
            if number_real_roots > 1
                x_ampl  = max(x_roots_real) - min(x_roots_real);
                x_LB    = min(x_roots_real) - 0.2*x_ampl;
                x_UB    = max(x_roots_real) + 0.2*x_ampl;
            elseif number_real_roots == 1
                x_UB    = x_roots_real + 1;
                x_LB    = x_roots_real - 1;
            else
                x_LB = -1;
                x_UB = 1;
            end
            
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    

            hold on
            grid on

            % Polynomial
            x_list  = linspace(x_LB, x_UB, 1e3);
            f_list  = a*x_list.^3 + b*x_list.^2 + c*x_list + d;
            plot(x_list, f_list, 'LineWidth', 2, 'color', 'b', 'DisplayName', 'Polynomial');

            % Zero-line
            plot([x_LB, x_UB], [0, 0], 'LineWidth', 2, 'color', 'k', 'HandleVisibility', 'Off');

            % Roots
            for k = 1 : number_real_roots
                x_root = x_roots_real(k);
                sc_root = scatter(x_root, 0, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Root');

                if k > 1
                    sc_root.HandleVisibility = 'Off';
                end
            end

            % Axes
            xlim([x_LB, x_UB]);
            xlabel('x');
            ylabel('f(x)');

            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off   
            
            % Pause message
            disp('The roots to the cubic polynomial have been found. The figure will close and script will finish upon a key-press');
            pause();
            
            close(1);
        end
end