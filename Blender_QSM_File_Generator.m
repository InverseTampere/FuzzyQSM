% This script takes a TreeQSM generated QSM stored in a .mat file and converts it into a .txt file for Markku Åkerblom's Blender add-on

function Blender_QSM_File_Generator(QSM, QSM_file_name)

%% Inputs %%
    % Cylinder data
    branch_index_list   = QSM.cylinder.branch;
    branch_index_list   = single(branch_index_list);                    % So that everything is in single format

    start_point_matrix  = QSM.cylinder.start;
    cyl_axis_matrix     = QSM.cylinder.axis;
    cyl_length_list     = QSM.cylinder.length;
    cyl_radius_list     = QSM.cylinder.radius;

%% Write the .txt file %%
    % Branch index - Start point - Axis - Length - Radius
    QSM_data_matrix = [branch_index_list, start_point_matrix, cyl_axis_matrix, cyl_length_list, cyl_radius_list];

    % Writing to the new .txt file
    QSM_file_name_parts = strsplit(QSM_file_name, '.');
    QSM_file_name_basis = QSM_file_name_parts{1};
    QSM_file_name       = [QSM_file_name_basis, '.txt'];

    writematrix(QSM_data_matrix, QSM_file_name, 'Delimiter', ' ');