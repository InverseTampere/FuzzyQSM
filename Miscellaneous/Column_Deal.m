% A deal-esque function that splits an input into its columns, 
% i.e. [x, y, z] = Column_Deal(vector) instead of [x, y, z] = deal(vector(1), vector(2), vector(3))

% If the input is a 1 x n cell array it returns the entries within the cell array and not individual cells

function varargout = Column_Deal(input_matrix)
    % Check whether the script can execute
    number_outputs                  = nargout;
    [number_rows, number_columns]   = size(input_matrix);

    if number_outputs > number_columns
        error('The number of expected outputs exceeds the number of columns of the input.');
    end

    % Assigning the outputs
    varargout = cell(1, number_outputs);

    for i = 1 : number_outputs
        % The i'th column is selected as the i'th output
        if iscell(input_matrix) && number_rows == 1
            input_column = input_matrix{i};             % It returns the content of the cell rather than a 1 x 1 cell array
        else
            input_column = input_matrix(:, i);          
        end

        % It is assigned as output
        varargout{i} = input_column;
    end        
end