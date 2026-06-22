% The given boolean is applied to the fields of the structure whenever possible

function Structure = Structure_Boolean(Structure, bool_list)

    %% Applying the boolean %%
        % Going through each field in the structure
        fields          = fieldnames(Structure);
        number_fields   = length(fields);

        for f = 1 : number_fields
            % This field's data
            field       = fields{f};
            field_data  = Structure.(field);

            % Check if it is a numeric or cell array
            if ~isnumeric(field_data) & ~iscell(field_data)
                continue
            end

            % Applying the boolean if possible
            [number_rows, number_columns] = size(field_data);

            if number_rows == numel(bool_list)
                field_data = field_data(bool_list, :);
            elseif number_columns == numel(bool_list)
                field_data = field_data(:, bool_list);
            end

            Structure.(field) = field_data;
        end
end