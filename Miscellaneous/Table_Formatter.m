% This script can format the given data to the desired format for the table, as well as write an .xls file if a file name is given (can be without .xls)
% When the table is displayed it shows entries as {'123.456'}, but the written file shows numbers correctly as 123.456

% Format options are the maximum number of digits and whether it is 'Exponential', 'Integer', or 'Float'
% Note that if the number exceeds the maximum, exponential notation is used regardless of other settings

function Table = Table_Formatter(data_array, number_digits, output_format, row_names, column_names, file_name, Print)

    %% Changing the format %%
        % Conversion to exponential format (also used if there are too many digits before the period)
        exp_formatting_string   = sprintf('%%.%ie', number_digits - 1);                             % Minus one, due to the digit in front        
        Exp_Formatting_fun      = @(data) sprintf(exp_formatting_string, data);
            
        exp_cell_array          = arrayfun(Exp_Formatting_fun, data_array, 'UniformOutput', false);

        % Conversion to integer format (used to count the number of digits before the period)
        Int_Formatting_fun      = @(data) num2str(round(data));
        int_cell_array          = arrayfun(Int_Formatting_fun, data_array, 'UniformOutput', false);
                
        % Checking whether there are too many digits before the period
        Number_Digits_fun       = @(integer_string) length(strrep(integer_string, '-', ''));        % The minus sign is not counted
        number_digits_bp_matrix = cellfun(Number_Digits_fun, int_cell_array);
        too_long_bool           = number_digits_bp_matrix > number_digits;

        % Conversion to exponential format
        if strcmp(output_format, 'Exponential')
            % The exponential entries are the desired output
            cell_array = exp_cell_array;

        % Conversion to integers
        elseif strcmp(output_format, 'Integer')
            % When they are too long, the exponential format is used instead
            cell_array                  = int_cell_array;
            cell_array(too_long_bool)   = exp_cell_array(too_long_bool);

        % Conversion to floats
        elseif strcmp(output_format, 'Float')
            % Function to convert to floats of varying precision

            Flt_Formatting_String_fun       = @(number_digits_bp) sprintf('%%.%if', number_digits - number_digits_bp);
            float_formatting_string_cell    = arrayfun(Flt_Formatting_String_fun, number_digits_bp_matrix, 'UniformOutput', false);

            Flt_Formatting_fun      = @(float_formatting_string, data) sprintf(float_formatting_string, data);
            data_cell_array         = num2cell(data_array);
            cell_array              = cellfun(Flt_Formatting_fun, float_formatting_string_cell, data_cell_array, 'UniformOutput', false);

            % When they are too long, the exponential format is used instead
            cell_array(too_long_bool) = exp_cell_array(too_long_bool);
        else
            error('The output format %s was not understood. It must be Exponential, Integer or Float', output_format);
        end

    %% Creating the table %%
        % The table is created from the array
        Table = array2table(cell_array, 'RowNames', row_names, 'VariableNames',  column_names);

        % Shown if desired
        if Print == true
            disp(Table);
        end

        % And written to a file if the file name is given
        if ~isempty(file_name)
            % The file format is given
            file_name = strrep(file_name, '.xls', '');          % It is removed if it was there already
            file_name = sprintf('%s.xls', file_name);           % And added

            writetable(Table, file_name, 'WriteRowNames', true);

            if Print == true
                fprintf('The table was saved as %s \n', file_name);
            end
        end
end