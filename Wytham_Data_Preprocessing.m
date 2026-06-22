% This script collects the original Wytham QSMs located in the folders in Data/Original_QSMs,
% converts them to the TreeQSM 2.4.1 structure which are saved in Data/Original_QSMs
% and creates .txt files for Blender which are saved in Data/Original_QSMs

clear variables
close all
clc

%% Inputs %%
    Point_Cloud_Plot        = true;             % [true, false] If desired, the point cloud can be checked to gauge the coverage

%% Conversion of QSMs %%
    % Original QSM directory files
    Original_Files  = dir('Data/Original_QSMs/**');
    number_files    = length(Original_Files);

    for f = 1 : number_files
        % File and folder names
        file_name   = Original_Files(f).name;
        folder_name = Original_Files(f).folder;
 
        if contains(file_name, 'wytham') && contains(file_name, '.mat')          % Check that it's a Wytham woods .mat file
            % Tree species and ID
            folder_parts    = strsplit(folder_name, '\');
            tree_species    = folder_parts{end};

            file_name_parts = strsplit(file_name, {'_', '-'});
            tree_ID         = file_name_parts{3};

            QSM_file_name   = sprintf('%s_%s_QSM', tree_species, tree_ID);

            % Retrieve and save the QSM
            QSM_file_loc    = sprintf('%s%s%s', folder_name, '\', file_name);
            QSM             = Wytham_QSM_Reader(QSM_file_loc, Point_Cloud_Plot);
            
            save(QSM_file_name, 'QSM');

            % Create a .txt file for Blender
            Blender_QSM_File_Generator(QSM, QSM_file_name);

            % Move the files to the correct folder
            movefile(sprintf('%s*', QSM_file_name), folder_name);
        end
    end
    