% The roots of the quartic polynomial a*x^4 + b*x^3 + c*x^2 + d*x + e are found using Ferrari's solution method
% They are split into a real and complex list, where the number of complex roots is always four

function [x_roots_real, number_real_roots, x_roots_complex] = Quartic_Polynomial(a, b, c, d, e, Plot, Print, Diagnostics)

    %% Roots %%
        % If a is nonzero, the 4th order polynomial is solved
        if a ~= 0
            % Converted to the depressed polynomial where a = 1
            B = b/a;
            C = c/a;
            D = d/a;
            E = e/a;
            
            % Change of variable to y = x - b/4 to remove cubic term (f(y) = y^4 + py^2 + qy + r = 0)
            p = (8*C - 3*B^2)/8;
            q = (B^3 - 4*B*C + 8*D)/8;
            r = (-3*B^4 + 256*E - 64*B*D + 16*B^2*C)/256;
            
            % If q = 0, f(y) is biquadratic and the roots can be found simply as z = y^2
            if abs(q) < 1e-12       % Matlab rounding can cause it to be slightly nonzero 
                z_roots             = 1/2 * (-p + [1, -1] * sqrt(p^2 - 4*r));
                x_roots_complex     = [1; -1] .* sqrt(z_roots) - B/4;
    
            % Otherwise, Ferrari's solution method is used
            else
                % The first real root to the cubic polynomial 8m^3 + 8pm^2 + (2p^2 - 8r)m - q^2 are found
                alpha   = 8;
                beta    = 8*p;
                gamma   = 2*p^2 - 8*r;
                delta   = -q^2;
                
                [m_roots_real, ~, ~] = Cubic_Polynomial(alpha, beta, gamma, delta, Diagnostics, Diagnostics);
                m = max(m_roots_real);
                
                % The roots are then
                x_roots_complex = 1/2 * ([1, -1]*sqrt(2*m) + [1; -1]*sqrt(-(2*p + 2*m + [1, -1]*sqrt(2/m)*q))) - B/4;
            end
            
            % The roots are put into column-format and real roots are selected
            x_roots_complex     = x_roots_complex(:);
            
            x_roots_real        = x_roots_complex(abs(imag(x_roots_complex)) < 1e-6);          % A very minor imaginary component may be present due to Matlab rounding
            x_roots_real        = real(x_roots_real);
            number_real_roots   = length(x_roots_real);
        
        % If a is zero, it is a cubic polynomial
        else
            [x_roots_real, number_real_roots, x_roots_complex] = Cubic_Polynomial(b, c, d, e, Diagnostics, Diagnostics);
        end

    %% Printed messages %%
        if Print == true
            fprintf('Roots: (%g are real) \n', number_real_roots);
            
            for r = 1 : 4
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
            f_list  = a*x_list.^4 + b*x_list.^3 + c*x_list.^2 + d*x_list + e;
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
            disp('The roots to the quartic polynomial have been found. The figure will close and script will finish upon a key-press');
            pause();
            
            close(1);
        end
end