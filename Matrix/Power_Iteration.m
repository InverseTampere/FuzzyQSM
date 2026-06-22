% This script performs power iteration to find the strongest and weakest eigenvectors
% The convergence threshold is expected to be a percentage difference between two iterations

function [strongest_eigenvector, weakest_eigenvector] = Power_Iteration(covariance_matrix, convergence_threshold)
    
    %% Inputs %%
        number_check_steps  = 1e1;      % The number of steps taken before checking for convergence
                                        % Usually 10 steps are sufficient for convergence, and checking after each step is inefficient
        max_steps           = 1e2;      % The maximum number of steps 

    %% Power iteration %%                                
        % Random vectors are used to initiate the iteration
        num_dim                 = size(covariance_matrix, 1);
        strongest_eigenvector   = rand(num_dim, 1);        
        weakest_eigenvector     = rand(num_dim, 1);
        
        % Inputs for the loop
        difference  = Inf;
        steps       = 0;
        
        while difference > convergence_threshold && steps <= max_steps
            steps = steps + number_check_steps;
            
            for i = 1 : number_check_steps
                strongest_eigenvector   = covariance_matrix * strongest_eigenvector;            % Extended by the covariance matrix (dispersion)
                strongest_eigenvector   = strongest_eigenvector / norm(strongest_eigenvector);
                
                weakest_eigenvector     = covariance_matrix \ weakest_eigenvector;              % Extended by the inverted covariance matrix (clustering)
                weakest_eigenvector     = weakest_eigenvector / norm(weakest_eigenvector);
            end
            
            % To compute the difference, the opposite matrix operations are performed to retrieve the previous vectors
            previous_strongest_eigenvector  = covariance_matrix \ strongest_eigenvector;
            previous_strongest_eigenvector  = previous_strongest_eigenvector / norm(previous_strongest_eigenvector);
            
            previous_weakest_eigenvector    = covariance_matrix * weakest_eigenvector;
            previous_weakest_eigenvector    = previous_weakest_eigenvector / norm(previous_weakest_eigenvector);
            
            difference_strong   = (strongest_eigenvector - previous_strongest_eigenvector) * 100;
            difference_weak     = (weakest_eigenvector - previous_weakest_eigenvector) * 100;
            
            % The greatest norm of the difference vector is used to assess convergence
            difference = max(norm(difference_strong), norm(difference_weak));
        end
        
        % Row-vectors are returned
        strongest_eigenvector   = strongest_eigenvector';
        weakest_eigenvector     = weakest_eigenvector';
        
end