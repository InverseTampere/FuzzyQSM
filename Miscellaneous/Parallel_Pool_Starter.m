% This script starts the parallel pool with the desired idle time until deactivation (timeout) and number of cores
% If the number of cores aren't specified, the maximum will be used

function Parallel_Pool_Starter(idle_timeout, number_cores)
    
    % The parallel pool is attempted to be started
    try
        if nargin == 0
            idle_timeout = 30;
            number_cores = feature('numcores');
        end
        if nargin == 1
            number_cores = feature('numcores');
        end
        
        parpool('local', number_cores, 'IdleTimeout', idle_timeout);
        
    % Warning if it could not be started and wasn't already
    catch
        if isempty(gcp('nocreate'))
            error('Parallel pool could not be started.')
        end
    end
end