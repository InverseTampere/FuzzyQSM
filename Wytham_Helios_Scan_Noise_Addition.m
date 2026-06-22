% This takes the Helios scans performed without noise and adds Gaussian noise to them, as well as removes points that exceed the maximum incidence angle

clear variables
close all
beep off
clc

%% Inputs %%
    % Scanner parameters
    beam_divergence             =   0.30;                  % [mrad] 2 sigma (1/e2) divergence half-angle of the beam
    beam_exit_diameter          =   0.00;                   % [mm] Diameter of +/-2 sigma of the beam at the scanner exit
    sigma_range_0               =   03.0;                   % [mm] The range uncertainty of the instrument itself
    range_bias                  =   00.0;                   % [mm] The range bias
    max_incidence_angle         =   80;                     % [deg] Max incidence angle can be limited to prevent points that are near-oblique, 0 to 90

    % File locations
    Helios_scans_folder         =   "Data/Helios_Scans";         % Folder in which the Helios scans are located
    smoothed_models_folder      =   "Data/Smoothed_Models";      % Folder in which the smoothed models are located
    output_folder               =   "Data/Point_Clouds/030_beamdivergence";         % Folder where the noisy point clouds will be stored

    % Diagnostic plots
    Plot                        =   false;                  % [true, false] Shows the original and noisy point clouds

%% Input structures %%
    % Conversion to SI units 
    beam_divergence             =   beam_divergence * 1e-3;            % mrad to rad
    beam_exit_diameter          =   beam_exit_diameter * 1e-3;         % mm to m
    sigma_range_0               =   sigma_range_0 * 1e-3;              % mm to m
    range_bias                  =   range_bias * 1e-3;                 % mm to m
    max_incidence_angle         =   deg2rad(max_incidence_angle);      % deg to rad

    % Scanner parameter structure
    Scanner_Parameters = struct('beam_divergence', beam_divergence, 'beam_exit_diameter', beam_exit_diameter, 'sigma_range_0', sigma_range_0, 'range_bias', range_bias, 'max_incidence_angle', max_incidence_angle);

    % Scanning parameter structure 
    scanning_parameter_file_path    = sprintf('%s%sScanning_Parameters.mat', Helios_scans_folder, '/');
    Scanning_Parameters             = load(scanning_parameter_file_path);

%% Adding noise to each tree's point cloud %%
    % Tree ID folders which contain the data
    tree_scan_folders   = dir(Helios_scans_folder);
    tree_scan_folders   = {tree_scan_folders(:).name};

    dot_folders         = contains(tree_scan_folders, '.');       % Remove the . and .. folders and Scanning_Parameters.mat
    tree_ID_cell        = tree_scan_folders(~dot_folders);

    number_trees = length(tree_ID_cell);

    for t = 1 : number_trees
        % The tree's ID
        tree_ID = tree_ID_cell{t};

        % Helios scan data
        tree_scan_folder    = sprintf('%s%s%s', Helios_scans_folder, '/', tree_ID);

        coarsening_factor   = 1;
        PC_Diagnostics      = false;
        Point_Cloud_Data    = Helios_Data_Reader(tree_scan_folder, coarsening_factor, PC_Diagnostics);
        point_cloud_cell    = Point_Cloud_Data.point_cloud_cell;

        % The smoothed model's mesh
        mesh_file_path      = sprintf('%s%s%s%s%s_Mesh.mat', smoothed_models_folder, '/', tree_ID, '/', tree_ID);
        Triangular_Mesh     = load(mesh_file_path);

        % Determine which triangle each point corresponds to
        triangle_index_file_name = sprintf('%s%s%s%s%s_Triangle_Indices.mat', smoothed_models_folder, '/', tree_ID, '/', tree_ID);

        if exist(triangle_index_file_name, 'file')
            Triangle_Index_File = load(triangle_index_file_name);
            triangle_index_cell = Triangle_Index_File.triangle_index_cell;
        else
            triangle_index_cell = Helios_Scan_Triangle_Association(point_cloud_cell, Scanning_Parameters, Triangular_Mesh);
            save(triangle_index_file_name, 'triangle_index_cell');
        end
        
        % Add noise
        Point_Cloud_Data_n = Helios_Scan_Noise_Addition(Point_Cloud_Data, Scanner_Parameters, Scanning_Parameters, Triangular_Mesh, triangle_index_cell, Plot);

        % Saving the result and including the scanning and scanner parameters
        point_cloud_file_name = sprintf('%s_Point_Cloud.mat', tree_ID);
        save(point_cloud_file_name, 'Point_Cloud_Data_n', 'Scanning_Parameters', 'Scanner_Parameters');
        movefile(point_cloud_file_name, output_folder);

        % Progress message
        fprintf('%i/%i: Noise has been added to the point cloud of %s. \n', t, number_trees, tree_ID);
    end