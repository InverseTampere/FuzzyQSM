% This script performs Helios scanning over all the tree objects

clear variables
close all
beep off
clc

%% Inputs %%
    %--% Scanning parameters %--%
    Scanner_loc_cell    = {[20, 0, 1.5], [0, 10, 1.5]};     % {[x,y,z], [x,y,z], ...} Locations of the scanners
    beam_divergence     = 0.60;                         % [mrad] Beam divergence full-angle, 1/e2

    %--% File locations %--%
    Helios_folder       = "C:\Users\svvive\AppData\Local\helios";
    Helios_batch_file   = "Menu\Helios++ Survey Session.bat";
    survey_file         = "data\surveys\Wytham_Woods\Wytham_Survey.xml";
    scene_file          = "data\scenes\Wytham_Woods\Wytham_Scene.xml";
    scanner_file        = "Lib\site-packages\pyhelios\data\scanners_tls.xml";
    scanner_ID          = "riegl_vz400_adj";
    tree_models_folder  = "C:\Users\svvive\OneDrive - TUNI.fi\Remote_Sensing_Forests\FuzzyQSM\Data\Smoothed_Models";
    output_folder       = "C:\Users\svvive\OneDrive - TUNI.fi\Remote_Sensing_Forests\FuzzyQSM\Data\Helios_Scans";

    %--% Outputs %--%
    Plot                = true;                        % [true, false]

%% Helios scanning %%
    % Tree models
    tree_model_folders  = dir(tree_models_folder);
    tree_model_folders  = {tree_model_folders(:).name};

    dot_folders         = contains(tree_model_folders, '.');       % Remove the . and .. folders and any files
    tree_model_folders  = tree_model_folders(~dot_folders);
    number_trees        = length(tree_model_folders);

    % Scanning parameter structure
    beam_divergence     = beam_divergence * 1e-3;
    Scanning_Parameters = struct('Scanner_loc_cell', {Scanner_loc_cell}, 'number_scanners', length(Scanner_loc_cell), 'beam_divergence', beam_divergence);

    scanning_parameters_file_name = sprintf('%s%s%s', output_folder, '\', 'Scanning_Parameters.mat');
    save(scanning_parameters_file_name, '-struct', 'Scanning_Parameters');

    % Scanning each tree
    for t = 1 : number_trees
        % This tree's .obj file
        tree_model_folder   = tree_model_folders{t};
        Tree_File           = dir(sprintf("%s%s%s%s*.obj", tree_models_folder, '\', tree_model_folder, '\'));
        tree_file           = sprintf("%s%s%s", Tree_File.folder, '\', Tree_File.name);

        tic;
        Running_Helios(Helios_folder, scene_file, survey_file, scanner_file, scanner_ID, tree_file, Helios_batch_file, output_folder, Scanning_Parameters, Plot);
        t_Helios = toc;

        current_time = datetime('now', 'format', 'dd_HHmmss');
        fprintf('t = %s. Helios scanning %s took %.3g s \n', current_time, Tree_File.name, t_Helios);
    end