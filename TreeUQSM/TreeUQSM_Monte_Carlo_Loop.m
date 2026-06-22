% TreeQSM is ran in a Monte Carlo loop by sampling the fuzzy cloud

function [MC_QSM_cell, TreeQSM_GMM] = TreeUQSM_Monte_Carlo_Loop(Point_Cloud_Distributions, Parallel_Pool, Fitting_Parameters, Statistical_Values, Monte_Carlo_Inputs, QSM_init, TreeQSM_Inputs, Output_Decisions)

    %% Structure inputs %%
        % Distributions
        distribution_mu_cell        = Point_Cloud_Distributions.distribution_mu_cell;

        % Parallel pool
        Parallel_Loop               = Parallel_Pool.Parallel_Loop;
        idle_timeout                = Parallel_Pool.idle_timeout;

        % Fitting parameters
        number_models               = Fitting_Parameters.number_models;
        Optimum_metric              = Fitting_Parameters.Optimum_metric;

        % Statistical values
        Bhatt_coeff_threshold       = Statistical_Values.Bhatt_coeff_threshold;
        
        % Monte Carlo simulation
        TreeQSM_Stochasticity       = Monte_Carlo_Inputs.TreeQSM_Stochasticity;
        Resampling                  = Monte_Carlo_Inputs.Resampling;
        max_MC_length               = Monte_Carlo_Inputs.max_MC_length;  
        number_GMM_fitment_iter     = Monte_Carlo_Inputs.number_GMM_fitment_iter;
        GMM_batch_size              = Monte_Carlo_Inputs.GMM_batch_size;

        % Initial QSM
        QSM_treedata                = QSM_init.treedata;

        % Outputs
        Print                       = Output_Decisions.Print;

    %% Manual inputs %%
        Sampling_Diagnostics        = false;        % [true, false]

    %% Preliminary %%
        % Number of loops
        number_loops    = floor(max_MC_length / GMM_batch_size);
        number_loops    = max(1, number_loops);
        MC_loop_length  = GMM_batch_size * number_loops;   

        % Initiate the parallel pool
        if Parallel_Loop == true
            % The parallel pool is started
            number_cores = feature('numcores');    
            Parallel_Pool_Starter(idle_timeout, number_cores)
        else
            % Otherwise the number of cores used in parallel are set to 0 to make the parfor run as a for
            MC_loop_length  = max_MC_length;
            number_cores    = 0;
        end

        % If TreeQSM stochasticity is not desired, the cover sets and segments are retrieved
        if TreeQSM_Stochasticity == false
            Cover_Sets  = Monte_Carlo_Inputs.Cover_Sets;
            Segments    = Monte_Carlo_Inputs.Segments;
        else
            Cover_Sets  = [];
            Segments    = [];
        end

    %% Monte Carlo loop %%
        % Each iteration's QSM is saved
        MC_QSM_cell = cell(1, MC_loop_length);

        % Gaussian mixture model
        tree_parameter_fields   = {'TrunkVolume', 'BranchVolume', 'BranchLength'};
        number_fields           = length(tree_parameter_fields);
        GMM_Data_Struct         = cell2struct(repmat({[]}, 1, number_fields), tree_parameter_fields, 2);
        TreeQSM_GMM             = struct('GMM', [], 'Data', GMM_Data_Struct);

        % Point cloud matrix
        point_cloud_matrix  = vertcat(distribution_mu_cell{:});

        % Progress counter
        DQ      = parallel.pool.DataQueue;
        tick    = 0;
        N       = MC_loop_length;
        afterEach(DQ, @ProgressUpdate);

        Convergence     = false;
        loop_counter    = 0;

        while Convergence == false && loop_counter < number_loops
            % The loop counter is updated until convergence or the maximum is reached
            loop_counter    = loop_counter + 1;
            
            start_ind       = (loop_counter - 1) * GMM_batch_size + 1;
            end_ind         = loop_counter * GMM_batch_size;

            %--% TreeQSM %--%      
            GMM_Variable_Data = TreeQSM_GMM.Data;               % Temporary structure which can be used in the parallel loop

            parfor (i = start_ind : end_ind, number_cores) 
                % Set the RNG seed for consistency (in a parallel loop the initially specified RNG seed is ignored)
                rng(i + 100);
    
                % Sampling the fuzzy cloud
                if Resampling == true
                    MC_Point_Cloud_Coord    = Gaussian_Point_Sampling(Point_Cloud_Distributions, Sampling_Diagnostics);
                    MC_point_cloud_cell     = MC_Point_Cloud_Coord.point_cloud_cell;
                    MC_point_cloud_matrix   = vertcat(MC_point_cloud_cell{:});
                else
                    MC_point_cloud_matrix   = point_cloud_matrix;
                end

                % Fitting the QSM
                if TreeQSM_Stochasticity == true
                    QSM_cell = cell(1, number_models);
    
                    for j = 1 : number_models
                        QSM_j       = treeqsm(MC_point_cloud_matrix, TreeQSM_Inputs);
                        QSM_cell{j} = QSM_j;
                    end
                    
                    MC_QSMs             = horzcat(QSM_cell{:});
                    [~, ~, ~, MC_QSM]   = select_optimum(MC_QSMs, Optimum_metric);
                else
                    MC_QSM = TreeQSM_Cylinder_Fitting(MC_point_cloud_matrix, Cover_Sets, Segments, TreeQSM_Inputs);
                end

                MC_QSM_cell{i} = MC_QSM;

                % The normalised Gaussian mixture model data
                QSM_Tree_Parameter_fun  = @(tree_parameter) MC_QSM.treedata.(tree_parameter) / QSM_treedata.(tree_parameter);               %#ok<PFBNS> 
                QSM_tree_parameter_cell = cellfun(QSM_Tree_Parameter_fun, tree_parameter_fields, 'UniformOutput', false);
                GMM_Variable_Data(i)    = cell2struct(QSM_tree_parameter_cell, tree_parameter_fields, 2);

                % Progress update
                send(DQ, i);
            end

            %--% Gaussian mixture model fitment %--%
            % The original structure is changed according to the temporary one
            TreeQSM_GMM.Data = GMM_Variable_Data;

            % The QSM data thus far
            [current_QSM_matrix, ~] = Structure_Data_Concatenation(TreeQSM_GMM.Data);

            % The optimal models are fitted based on the AIC
            [GM_Model, Shared_Covariance, number_GM_components, AICc_min] = Gaussian_Mixture_Model_Fitting(current_QSM_matrix, number_GMM_fitment_iter, []);
            TreeQSM_GMM.GMM{loop_counter} = GM_Model;

            %--% Convergence check %--%
            if loop_counter > 1 && ~isempty(GM_Model)           % If a GMM could not be fitted, empty values are given
                % The previous loop's Gaussian mixture model
                GMM_prev = TreeQSM_GMM.GMM{loop_counter - 1};

                % The Bhattacharyya coefficient is computed between the current and previous models
                if ~isempty(GMM_prev)
                    [Bhatt_coeff, Bhatt_distance] = Bhattacharyya_GMM(GM_Model, GMM_prev);

                % If a previous model does not exist, NaN values are given
                else
                    [Bhatt_coeff, Bhatt_distance] = deal(NaN);
                end
            else
                [Bhatt_coeff, Bhatt_distance] = deal(NaN);
            end

            % All the GMM metrics are saved
            TreeQSM_GMM.GMM_Metrics(loop_counter) = struct('Shared_Covariance', Shared_Covariance, 'number_GM_components', number_GM_components, 'AICc_min', AICc_min, 'Bhatt_coeff', Bhatt_coeff, 'Bhatt_distance', Bhatt_distance);
                                          
            %--% Convergence check %--%            
            % It has converged if the coefficient is above the threshold            
            number_MC_iterations = end_ind;

            if Bhatt_coeff > Bhatt_coeff_threshold
                % Convergence is set to true and the remaining entries are removed
                Convergence = true;
                MC_QSM_cell(number_MC_iterations + 1 : end) = [];
                
                % Convergence information is printed if desired
                if Print == true
                    fprintf('\n');
                    disp('---------------------------');
                    fprintf('%g Monte Carlo iterations were performed \n', number_MC_iterations);
                    fprintf('The GMM has a Bhattacharyya coefficient of %.3g \n', Bhatt_coeff);
                    disp('---------------------------');
                    fprintf('\n');
                end
            end
        end

        % A warning or error is displayed if the model has not converged within the number of steps
        if Convergence == false
            if isnan(Bhatt_coeff)
                warning('The Gaussian mixture model could not be created within %g iterations. \n', max_MC_length);
            else
                warning('The algorithm could not converge within %g iterations. The Bhattacharyya coefficient is still %.3g \n', max_MC_length, Bhatt_coeff);
            end            
        end

    %% Local functions %%
        % Progress update
        function ProgressUpdate(~)
            % Ensures that at most every P percent is printed
            P = 1;

            % The tick is updated
            tick = tick + 1;    

            % The last tick's and current tick's progress relative to P
            progress_last   = (tick - 1) / N * 100 / P;
            progress        = tick / N * 100 / P;

            if floor(progress) - floor(progress_last) >= 1
                current_time = datetime('now', 'format', 'dd_HHmmss');
                fprintf('   t = %s. Monte Carlo QSM fitting progress: %i QSMs \n', current_time, tick);
            end            
        end

end