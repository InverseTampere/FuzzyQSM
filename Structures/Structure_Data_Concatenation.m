% The arrays of the given structure are concatenated column-wise. They are thus expected to have an identical number of rows.
% The column indices in the matrix of each field also returned.

function [data_matrix, Field_Indices] = Structure_Data_Concatenation(Data_Structure)

    %% Concatenated data %%
        % The data is first converted into a cell array and then into a matrix
        data_cell   = squeeze(struct2cell(Data_Structure))';
        data_matrix = cell2mat(data_cell);

    %% Indices of each field %%
        % The number of columns of each field
        Number_Columns_fun  = @(matrix) size(matrix, 2);
        number_columns_list = structfun(Number_Columns_fun, Data_Structure(1));             % The first entry is taken in case the structure has more than one entry

        % The indices of each field
        cumsum_columns_list = cumsum(number_columns_list);
        start_ind_list      = [1; cumsum_columns_list(1 : end - 1) + 1];

        Indices_fun         = @(number_columns, start_ind) start_ind : start_ind + number_columns - 1;
        field_indices_cell  = arrayfun(Indices_fun, number_columns_list, start_ind_list, 'UniformOutput', false);

        % Create a structure
        field_names     = fieldnames(Data_Structure);
        Field_Indices   = cell2struct(field_indices_cell, field_names, 1);
       
end