% This script reads the .obj files (which are relatively slow to read) and saves the data in .mat files

clear variables
close all
clc

%% Inputs %%
    tree_models_folder  = "C:\Users\svvive\OneDrive - TUNI.fi\Remote_Sensing_Forests\FuzzyQSM\Data\Smoothed_Models";        % Folder in which the smoothed models are located
    Obj_Diagnostics     = false;

%% File conversion %% 
    % Folders in which the data is located
    tree_model_folders  = dir(tree_models_folder);
    tree_model_folders  = {tree_model_folders(:).name};

    dot_folders         = contains(tree_model_folders, '.');       % Remove the . and .. folders and files
    tree_model_folders  = tree_model_folders(~dot_folders);
    number_trees        = length(tree_model_folders);

    for t = 1 : number_trees
        % This tree's .obj file
        tree_model_folder   = tree_model_folders{t};
        Tree_File           = dir(sprintf("%s%s%s%s*.obj", tree_models_folder, '\', tree_model_folder, '\'));
        folder_name         = Tree_File.folder;
        obj_file_name       = Tree_File.name;
        tree_obj_file_path  = sprintf("%s%s%s", Tree_File.folder, '\', Tree_File.name);

        % The vertices and mesh triangles in the data
        delimiter               = ' ';
        Y_Axis_Up               = true;                             
        [Vertices, Triangles]   = Object_File_Reader(tree_obj_file_path, delimiter, Y_Axis_Up, Obj_Diagnostics);

        % Save the results
        tree_file_name = sprintf('%s_Mesh.mat', tree_model_folder);
        save(tree_file_name, 'Vertices', 'Triangles');

        movefile(tree_file_name, folder_name);
    end