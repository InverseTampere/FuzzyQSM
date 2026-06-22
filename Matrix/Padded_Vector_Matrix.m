% This script takes unequal length vectors in a cell array, and combines them into one padded matrix 
% Each vector is given a separate column, and empty entries are denoted with the given padding value

function padded_matrix = Padded_Vector_Matrix(vector_cell, padding_value)
    
    %% The size of the padded matrix %%
        % The size is dictated by the number of vectors and their maximum length
        num_entries_list    = cellfun(@length, vector_cell);
        max_length          = max(num_entries_list);

    %% Filling the padded matrix %%
        % For cellfun, a cell array is created containing the number of entries        
        num_entries_cell    = num2cell(num_entries_list);
    
        % Vectors are returned in a cell array using cellfun with the following function
        Padding_fun         = @(vector, num_entries) [vector; repmat(padding_value, max_length - num_entries, 1)];
        padded_cell_array   = cellfun(Padding_fun, vector_cell, num_entries_cell, 'UniformOutput', false);
        
        % And then converted to a matrix
        padded_matrix = horzcat(padded_cell_array{:});

end