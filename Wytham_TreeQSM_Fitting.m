% This script fits QSMs to all the noisy point clouds using TreeQSM

clear variables
close all
beep off
clc

%% Inputs %%
    % Set rng
    rng(1);

    %--% Data locations %--%
    point_clouds_folder         =   'Data/Point_Clouds/030_beamdivergence';    % Folder in which the noisy point clouds are located
    discrete_QSMs_folder        =   'Data/TreeQSMs/030_beamdivergence';        % Folder in which the TreeQSM generated QSMs will be put

    %--% Data parameters %--%
    coarsening_factor           =   001;                    % [-] Integer by which the number of points is reduced

    %--% TreeQSM parameters %--%
    % Patch generation parameters
    PG_Inputs_Defined           =   true;                   % [true, false] If true, the upcoming parameter values are instead determined by define_inputs.m
    patch_diameter_1            =   [0.08, 0.12];           % [-] Patch size of the first uniform-size cover. Note that a vector of options may be given
    patch_diameter_2_min        =   [0.02, 0.03];           % [-] Minimum patch size for the second cover. Note that a vector of options may be given
    patch_diameter_2_max        =   [0.07, 0.10];           % [-] Maximum patch size for the second cover. Note that a vector of options may be given
    
    ball_radius_1_increase      =   0.015;                  % [-] The ball radius of the first cover is slightly larger than the first patch diameter
    min_number_ball_points_1    =   3;                      % [-] Minimum number of points for the first balls. Generally, 3 is a good value
    ball_radius_2_increase      =   0.010;                  % [-] The ball radius of the second cover is slightly larger than the maximum of the second patch diameter
    min_number_ball_points_2    =   1;                      % [-] Minimum number of points for the second balls. Generally, 1 is a good value

    % Cylinder radius correction
    min_cylinder_radius         =   0.0025;                 % [m] Minimum cylinder radius
    Parent_correction           =   true;                   % [true, false] The radius of a child should be smaller than its parent
    Taper_correction            =   true;                   % [true, false] Taper correction uses partially linear (stem) and parabolic (branch) taper correction

    % Growth-volume correction
    Growth_volume_correction    =   false;                  % [true, false] Whether or not the growth-volume correction approach by Jan Hackenberg is used
                                                            %               If true, parent and taper correction are suggested to be off
    growth_volume_factor        =   2.5;                    % [-] Defines upper and lower bounds

    % Optimal model selection
    number_models               =   05;                         % [-] The number of models that are fit each time to counter-act stochasticity
    Optimum_metric              =   'trunk+branch_mean_dis';    % Metric for selecting the optimal model. See select_optimum.m for the different options

    %--% Outputs %--%
    TreeQSM_Print               =   false;                  % [true, false] Prints TreeQSM's results
    TreeQSM_Plot                =   false;                  % [true, false] Plots TreeQSM's results

%% Structures %%
    % Fitting parameters for TreeQSM
    Fitting_Parameters  = struct('patch_diameter_1', patch_diameter_1, 'patch_diameter_2_min', patch_diameter_2_min, 'patch_diameter_2_max', patch_diameter_2_max, ...
                                 'ball_radius_1_increase', ball_radius_1_increase, 'min_number_ball_points_1', min_number_ball_points_1, 'ball_radius_2_increase', ball_radius_2_increase, 'min_number_ball_points_2', min_number_ball_points_2, ...
                                 'min_cylinder_radius', min_cylinder_radius, 'Parent_correction', Parent_correction, 'Taper_correction', Taper_correction, 'Growth_volume_correction', Growth_volume_correction, 'growth_volume_factor', growth_volume_factor, 'number_models', number_models, 'Optimum_metric', Optimum_metric);

    % Data parameters
    Data_Parameters     = struct('coarsening_factor', coarsening_factor, 'Only_tree', true);

%% Fitting the QSMs %%
    % Folders which contain the data
    point_cloud_files   = dir(point_clouds_folder);
    point_cloud_files   = {point_cloud_files(:).name};

    dot_folders         = startsWith(point_cloud_files, '.');       % Remove the . and .. folders
    point_cloud_files   = point_cloud_files(~dot_folders);

    number_trees = length(point_cloud_files);

    for t = 1 : number_trees
        % Read the data
        point_cloud_file    = point_cloud_files{t};
        Point_Cloud_File    = load(sprintf('%s/%s', point_clouds_folder, point_cloud_file));

        point_cloud_cell    = Point_Cloud_File.Point_Cloud_Data_n.point_cloud_cell;
        point_cloud_matrix  = vertcat(point_cloud_cell{:});

        % An input structure is created based on the inputs specified above
        tree_ID             = strrep(point_cloud_file, '_Point_Cloud.mat', '');
        TreeQSM_Inputs      = TreeQSM_Input_Creator(Fitting_Parameters, Data_Parameters, PG_Inputs_Defined, point_cloud_matrix, tree_ID, TreeQSM_Print, TreeQSM_Plot);

        tic;
        % The optimal discrete QSM is selected
        QSMs_discrete_cell = cell(1, number_models);
    
        for i = 1 : number_models
            QSMs_discrete_i         = treeqsm(point_cloud_matrix, TreeQSM_Inputs);
            QSMs_discrete_cell{i}   = QSMs_discrete_i;
        end            
        
        QSMs_discrete                       = horzcat(QSMs_discrete_cell{:});
        [~, ~, Opt_Inputs, QSM_Discrete]    = select_optimum(QSMs_discrete, Optimum_metric);
    
        % Non-existant branches are removed
        QSM_Discrete = TreeQSM_Branch_Trimmer(QSM_Discrete);

        % Parentless branches may erroneously exist, and are removed here
        [QSM_Discrete, ~] = QSM_Parentless_Branches(QSM_Discrete);                  % Note that the tree volumes are now in cubic metres, not litres

        % The QSM is saved
        discrete_QSM_file_name = sprintf('%s/%s_Discrete_QSM.mat', discrete_QSMs_folder, tree_ID);
        save(discrete_QSM_file_name, 'QSM_Discrete', 'QSMs_discrete_cell', 'Opt_Inputs', 'TreeQSM_Inputs', 'coarsening_factor');
        
        t_QSM = toc;
    
        current_time = datetime('now', 'format', 'dd_HHmmss');
        fprintf('t = %s. Discrete fitting %s took %.3g s \n', current_time, tree_ID, t_QSM);
    end