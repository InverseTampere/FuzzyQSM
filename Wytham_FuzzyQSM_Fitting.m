% This script fits QSMs to all the Helios scan data using a fuzzy approach
% Note that it expects the discrete QSMs to have already been fitted using TreeQSM

clear variables
close all
beep off
clc

%% Inputs %%
    %--% Execution %--% 
    rng(1);                                                 % Set the rng seed
    Parallel_Loop               =   false;                   % [true, false] Determines whether or not the Monte Carlo loop is ran in parallel
    max_cores                   =   10;                     % [-] If empty, the computer's number of cores is used. Otherwise this integer
    idle_timeout                =   030;                    % [min] Time until the parallel pool shuts down if idle        

    %--% Data %--%
    coarsening_factor           =   001;                    % [-] Integer by which the number of points is reduced
    discrete_QSMs_folder        =   'Data/TreeQSMs/030_beamdivergence';        % String of the folder in which the discrete QSMs are saved
    point_clouds_folder         =   'Data/Point_Clouds/030_beamdivergence';    % Folder in which the noisy point clouds are located
    fuzzy_QSMs_folder           =   'Data/FuzzyQSMs/030_beamdivergence';       % String of the folder in which the fuzzy QSMs are saved

    %--% Fuzzy fitting parameters %--%
    Fuzzy_Vector                =   true;                  % [true, false] Whether or not the vector is fitted or just taken from the discrete fit
    bounds_margin               =   0.50;                   % [-] Factor by which the geometry parameters can deviate from the initial estimates

    %--% Outputs %--%
    Print                       =   true;                  % [true, false] Shows printed statements regarding intermediate and final results
    Plot                        =   true;                  % [true, false] Shows plots of final results
    Compute_Radius_Errors       =   true;                  % [true, false] Computing the radius errors is slow, so can be toggled off

%% Input structures %%
    %--% Creation of the structures %--%
    Parallel_Pool           = struct('Parallel_Loop', Parallel_Loop, 'max_cores', max_cores, 'idle_timeout', idle_timeout);
    Fitting_Parameters      = struct('Fuzzy_Vector', Fuzzy_Vector, 'bounds_margin', bounds_margin);

%% Fuzzy QSM fitting %%
    % Folders which contain the data
    point_cloud_files   = dir(point_clouds_folder);
    point_cloud_files   = {point_cloud_files(:).name};

    dot_folders         = startsWith(point_cloud_files, '.');       % Remove the . and .. folders
    point_cloud_files   = point_cloud_files(~dot_folders);

    number_trees = length(point_cloud_files);

    for t = 1 : number_trees
        % Point cloud data
        point_cloud_file    = point_cloud_files{t};
        Point_Cloud_File    = load(sprintf('%s/%s', point_clouds_folder, point_cloud_file));
        Point_Cloud_Data    = Point_Cloud_File.Point_Cloud_Data_n;
        Scanner_Parameters  = Point_Cloud_File.Scanner_Parameters;
        Scanning_Parameters = Point_Cloud_File.Scanning_Parameters;

        Scanner_Parameters.sigma_range_device = Scanner_Parameters.sigma_range_0;

        % Discrete QSM fitted by TreeQSM
        tree_ID                 = strrep(point_cloud_file, '_Point_Cloud.mat', '');
        discrete_QSM_file_name  = sprintf('%s/%s_Discrete_QSM.mat', discrete_QSMs_folder, tree_ID);
        Discrete_QSM_File       = load(discrete_QSM_file_name);
        
        [QSM_Discrete, TreeQSM_Inputs, coarsening_factor_TreeQSM] = deal(Discrete_QSM_File.QSM_Discrete, Discrete_QSM_File.TreeQSM_Inputs, Discrete_QSM_File.coarsening_factor);

        if coarsening_factor_TreeQSM ~= coarsening_factor
            error('TreeQSM was used with a different coarsening factor so the QSMs are incompatible');
        end

        % Effectively the discrete cylinders are replaced by fuzzy cylinders
        tic;
        [QSM_Fuzzy, Cyl_Point_Cloud_Distributions_cell] = Fuzzy_QSM_Fitting(QSM_Discrete, Point_Cloud_Data, Fitting_Parameters, Scanner_Parameters, Scanning_Parameters, TreeQSM_Inputs, Parallel_Pool);
        t_QSM = toc;
    
        current_time = datetime('now', 'format', 'dd_HHmmss');
        fprintf('t = %s. Fuzzy fitting %s took %.3g s \n', current_time, tree_ID, t_QSM);
    
        % The distributions and fuzzy QSM are saved
        distributions_file_name = sprintf('%s/%s_Point_Cloud_Distributions.mat', fuzzy_QSMs_folder, tree_ID);
        save(distributions_file_name, 'Cyl_Point_Cloud_Distributions_cell');

        fuzzy_QSM_file_name = sprintf('%s/%s_Fuzzy_QSM.mat', fuzzy_QSMs_folder, tree_ID);
        save(fuzzy_QSM_file_name, 'QSM_Fuzzy');
    end