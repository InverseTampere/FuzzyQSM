% This script takes an initial QSM and fits fuzzy cylinders to each cylinder's point cloud

function [Fuzzy_QSM, Cyl_Point_Cloud_Distributions_cell] = Fuzzy_QSM_Fitting(QSM_init, Point_Cloud_Data, Fitting_Parameters, Scanner_Parameters, Scanning_Parameters, TreeQSM_Inputs, Parallel_Pool)

    %% Structure inputs %%
        % Initial QSM cylinders
        cylinder_indices_cell       = QSM_init.pmdistance.CylInd;
        Initial_Cylinders           = QSM_init.cylinder;
        init_cylinder_start_matrix  = Initial_Cylinders.start;
        init_cylinder_axis_matrix   = Initial_Cylinders.axis;
        init_cylinder_radius_list   = Initial_Cylinders.radius;
        init_cylinder_length_list   = Initial_Cylinders.length;
        init_cyl_branch_list        = Initial_Cylinders.branch;
        init_cyl_branch_order_list  = Initial_Cylinders.BranchOrder;
        init_cyl_branch_pos_list    = Initial_Cylinders.PositionInBranch;
        init_cyl_parent_list        = Initial_Cylinders.parent;
        init_cyl_surf_cov_list      = Initial_Cylinders.SurfCov;

        % Fitting parameters
        Fuzzy_Vector                = Fitting_Parameters.Fuzzy_Vector;

        % Scanning parameters
        Scanner_loc_cell            = Scanning_Parameters.Scanner_loc_cell;
        number_scanners             = Scanning_Parameters.number_scanners;

        % Point cloud
        point_cloud_cell            = Point_Cloud_Data.point_cloud_cell;
        number_points_list          = Point_Cloud_Data.number_points_list;

        % Parallel pool
        Parallel_Loop               = Parallel_Pool.Parallel_Loop;
        max_cores                   = Parallel_Pool.max_cores;
        idle_timeout                = Parallel_Pool.idle_timeout;

    %% Preliminary %%
        % Initiate the parallel pool
        if Parallel_Loop == true
            % The parallel pool is started
            if isempty(max_cores)
                number_cores = feature('numcores');
            else
                number_cores = max_cores;
            end
            
            Parallel_Pool_Starter(idle_timeout, number_cores)
        else
            % Otherwise the number of cores used in parallel are set to 0 to make the parfor run as a for
            number_cores = 0;
        end

        % Uniform point weights are used
        number_points                           = sum(number_points_list);
        uniform_weights_list                    = ones(number_points, 1);
        Fitting_Parameters.point_weights_list   = uniform_weights_list;

    %% Infinite cylinder fitting %%        
        % Each cylinder's point cloud
        point_cloud_matrix          = vertcat(point_cloud_cell{:});
        Cylinder_Point_Cloud_fun    = @(cylinder_indices) point_cloud_matrix(cylinder_indices, :);
        cylinder_point_cloud_cell   = cellfun(Cylinder_Point_Cloud_fun, cylinder_indices_cell, 'UniformOutput', false);

        cum_number_points_list      = [0, cumsum(number_points_list)];

        % The distributions are saved
        number_cylinders                    = length(cylinder_indices_cell);
        Cyl_Point_Cloud_Distributions_cell  = cell(1, number_cylinders);

        % Geometry parameters are set to the initial value by default, fuzziness is marked by the boolean
        fuzzy_cyl_bool_list     = false(number_cylinders, 1);

        init_cyl_centre_matrix  = init_cylinder_start_matrix + init_cylinder_length_list/2 .* init_cylinder_axis_matrix;
        fuzzy_cyl_centre_matrix = init_cyl_centre_matrix;
        fuzzy_cyl_radius_list   = init_cylinder_radius_list;
        fuzzy_cyl_axis_matrix   = init_cylinder_axis_matrix;
        fuzzy_cyl_length_list   = init_cylinder_length_list;

        % Progress counter
        DQ      = parallel.pool.DataQueue;
        tick    = 0;
        N       = number_cylinders;
        afterEach(DQ, @ProgressUpdate);

        parfor (c = 1 : number_cylinders, number_cores)     
            %--% Cylinder's data %--%
            % Geometry
            init_cylinder_start     = init_cylinder_start_matrix(c, :);
            init_cylinder_axis      = init_cylinder_axis_matrix(c, :);
            init_cylinder_radius    = init_cylinder_radius_list(c);
            init_cylinder_length    = init_cylinder_length_list(c);

            Init_Cylinder = struct('start', init_cylinder_start, 'axis', init_cylinder_axis, 'radius', init_cylinder_radius, 'length', init_cylinder_length);

            % Point cloud
            cylinder_point_cloud_total  = cylinder_point_cloud_cell{c};
            cylinder_indices_list       = cylinder_indices_cell{c};

            if isempty(cylinder_point_cloud_total)
                continue
            end

            % Outlier removal
            layer_height        = 0.02;         % Same as in cylinders.m in TreeQSM
            number_sectors      = 20;           % Same as in cylinders.m in TreeQSM
            [filter_bool, ~]    = surface_coverage_filtering(cylinder_point_cloud_total, Init_Cylinder, layer_height, number_sectors);

            cylinder_point_cloud_f          = cylinder_point_cloud_total(filter_bool, :);
            cylinder_point_cloud_cell{c}    = cylinder_point_cloud_f;

            % Point cloud per scanner
            cyl_point_cloud_cell    = cell(1, number_scanners);
            cyl_number_points_list  = zeros(1, number_scanners);

            for s = 1 : number_scanners
                ind_start   = cum_number_points_list(s) + 1;                %#ok<PFBNS>
                ind_end     = cum_number_points_list(s + 1);

                scanner_cyl_bool            = cylinder_indices_list >= ind_start & cylinder_indices_list <= ind_end;
                total_bool                  = scanner_cyl_bool & filter_bool;
                cyl_point_cloud_cell{s}     = cylinder_point_cloud_total(total_bool, :);
                cyl_number_points_list(s)   = sum(scanner_cyl_bool);
            end

            % Number of points per scanner's point cloud
            Number_Points_fun       = @(point_cloud) size(point_cloud, 1);
            cyl_number_points_list  = cellfun(Number_Points_fun, cyl_point_cloud_cell);

            if max(cyl_number_points_list) == 0
                continue
            end

            % Empty point clouds are removed
            empty_bool              = cyl_number_points_list == 0;
            cyl_number_points_list  = cyl_number_points_list(~empty_bool);
            cyl_point_cloud_cell    = cyl_point_cloud_cell(~empty_bool);

            Cyl_Point_Cloud_Coord   = struct('point_cloud_cell', {cyl_point_cloud_cell}, 'number_points_list', cyl_number_points_list);

            % As are associated scanners
            cyl_scanner_loc_cell    = Scanner_loc_cell(~empty_bool);        %#ok<PFBNS>
            cyl_number_scanners     = length(cyl_scanner_loc_cell);

            Cyl_Scanning_Parameters = struct('Scanner_loc_cell', {cyl_scanner_loc_cell}, 'number_scanners', cyl_number_scanners);

            %--% Fuzzy fitting %--%
            % Infinite cylinder fitting
            if Fuzzy_Vector == true
                % Fuzzy infinite cylinder
                [Fuzzy_Inf_Cylinder, Cyl_Point_Cloud_Distributions] = Fuzzy_QSM_Infinite_Cylinder_Fitting(Init_Cylinder, Cyl_Point_Cloud_Coord, Fitting_Parameters, Scanner_Parameters, Cyl_Scanning_Parameters);
            else
                % Fuzzy cross-section using the initial vector estimate
                [Fuzzy_Inf_Cylinder, Cyl_Point_Cloud_Distributions] = Fuzzy_QSM_Circle_Fitting(Init_Cylinder, Cyl_Point_Cloud_Coord, Fitting_Parameters, Scanner_Parameters, Cyl_Scanning_Parameters);
            end

            cylinder_centre = Fuzzy_Inf_Cylinder.centre;
            cylinder_radius = Fuzzy_Inf_Cylinder.radius;
            cylinder_axis   = Fuzzy_Inf_Cylinder.axis;

            % Assigning this cylinder's distributions
            Cyl_Point_Cloud_Distributions_cell{c} = Cyl_Point_Cloud_Distributions;

            % Cylinder top and bottom
            distribution_mu_matrix = vertcat(Cyl_Point_Cloud_Distributions.distribution_mu_cell{:});
            [projected_point_matrix, delta_list, ~] = Point_to_Vector_Projection(distribution_mu_matrix, cylinder_axis, cylinder_centre);
            
            [delta_bot, bot_ind]    = min(delta_list);
            cylinder_bottom         = projected_point_matrix(bot_ind, :);

            [delta_top, top_ind]    = max(delta_list);
            cylinder_top            = projected_point_matrix(top_ind, :);

            cylinder_length = delta_top - delta_bot;
            cylinder_centre = (cylinder_bottom + cylinder_top) / 2;

            % Inserting values into the arrays
            fuzzy_cyl_bool_list(c)          = true;
            fuzzy_cyl_radius_list(c)        = cylinder_radius;
            fuzzy_cyl_centre_matrix(c, :)   = cylinder_centre;
            fuzzy_cyl_axis_matrix(c, :)     = cylinder_axis;
            fuzzy_cyl_length_list(c)        = cylinder_length;

            % Progress update
            send(DQ, c);
        end

        % Cylinder structure
        Fuzzy_Cylinders = struct('radius', fuzzy_cyl_radius_list, 'length', fuzzy_cyl_length_list, 'centre', fuzzy_cyl_centre_matrix, 'axis', fuzzy_cyl_axis_matrix, 'branch', init_cyl_branch_list, 'BranchOrder', init_cyl_branch_order_list, 'PositionInBranch', init_cyl_branch_pos_list, 'SurfCov', init_cyl_surf_cov_list, 'parent', init_cyl_parent_list, 'fuzzy_bool', fuzzy_cyl_bool_list);      

        % Ensuring the branch cylinder axes are reasonable and point down the branch
        Fuzzy_Cylinders = Branch_Axis_Correction(Fuzzy_Cylinders, Initial_Cylinders);

        % Correct the radii
        Fuzzy_Cylinders = Branch_Radius_Correction(Fuzzy_Cylinders, TreeQSM_Inputs, cylinder_point_cloud_cell);

        % Ensuring the cylinder interfaces are reasonable
        Fuzzy_Cylinders = Branch_Interface_Fitting(Fuzzy_Cylinders);

    %% QSM structure %%
        % Cylinder volume
        fuzzy_cyl_length_list   = Fuzzy_Cylinders.length;
        fuzzy_cyl_radius_list   = Fuzzy_Cylinders.radius;
        fuzzy_cyl_volume_list   = pi*fuzzy_cyl_radius_list.^2 .* fuzzy_cyl_length_list;
        Fuzzy_Cylinders.volume  = fuzzy_cyl_volume_list;

        % Branch metrics
        [Branch_Metrics, Branch_Order_Metrics] = QSM_Branch_Metrics(Fuzzy_Cylinders, QSM_init);

        % Tree metrics
        Tree_Metrics = QSM_Tree_Metrics(Fuzzy_Cylinders, Branch_Order_Metrics);

        % QSM structure
        Fuzzy_QSM = struct('cylinder', Fuzzy_Cylinders, 'branch', Branch_Metrics, 'branch_order', Branch_Order_Metrics, 'treedata', Tree_Metrics);

    %% Local functions %%
        % Progress update
        function ProgressUpdate(~)
            % Ensures that at most every P percent is printed
            P = 1;

            % The tick is updated
            tick = tick + 1;    

            % The last tick's and current tick's progress relative to P
            progress_last   = floor((tick - 1) / N * 100 / P);
            progress        = floor(tick / N * 100 / P);

            if progress - progress_last >= 1
                current_time = datetime('now', 'format', 'dd_HHmmss');
                fprintf('   t = %s. Monte Carlo cylinder fitting progress: %i %% \n', current_time, progress);
            end            
        end

end