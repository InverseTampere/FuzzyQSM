% The tree metrics are determined of the Wytham smoothed tree models

clear variables
close all
beep off
clc

%% Inputs %%
    % Data locations
    smoothed_models_folder  = "C:\Users\svvive\OneDrive - TUNI.fi\Remote_Sensing_Forests\FuzzyQSM\Data\Smoothed_Models";
    original_QSMs_folder    = "C:\Users\svvive\OneDrive - TUNI.fi\Remote_Sensing_Forests\FuzzyQSM\Data\Original_QSMs";

    % Execution
    Parallel_Loop           = true;     % [true, false]
    max_cores               = 10;       % [-]
    idle_timeout            = 30;       % [min]
        
%% Tree metrics %%
    % Structure for parallel execution
    Parallel_Pool = struct('Parallel_Loop', Parallel_Loop, 'max_cores', max_cores, 'idle_timeout', idle_timeout);

    % Tree ID folders which contain the data
    smoothed_model_folders  = dir(smoothed_models_folder);
    smoothed_model_folders  = {smoothed_model_folders(:).name};

    dot_folders             = contains(smoothed_model_folders, '.');       % Remove the . and .. folders
    tree_ID_cell            = smoothed_model_folders(~dot_folders);

    number_trees = length(tree_ID_cell);

    for t = 1 : number_trees
        % This tree's data
        tree_ID                 = tree_ID_cell{t};

        QSM_file                = sprintf('%s/%s/%s_QSM.mat', original_QSMs_folder, tree_ID, tree_ID);
        smoothed_model_obj_file = sprintf('%s/%s/%s.obj', smoothed_models_folder, tree_ID, tree_ID);

        % Its metrics
        [Smoothed_Model_Metrics, Smoothed_Model_Circle_Geometry, Original_QSM] = Smoothed_Model_Tree_Metrics(QSM_file, smoothed_model_obj_file, Parallel_Pool);

        metrics_file_name = sprintf('%s/%s/%s_Metrics.mat', smoothed_models_folder, tree_ID, tree_ID);

        fprintf('The metrics have been determined for %i trees \n', t);
    end