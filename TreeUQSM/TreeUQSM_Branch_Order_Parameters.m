% This script aggregates the mean values and sums for the different branch parameters and returns
% these lists as well as their mean and standard deviations in structure arrays of form [scanner_combinations, branch_orders]

% Either the initial or MC QSM cell may be empty

function [Branch_Order_Parameters, Number_Branches, max_branch_order] = TreeUQSM_Branch_Order_Parameters(QSM_init_cell, MC_QSM_cell, branch_parameters)

    %% Maximum branch order %%
        % Retrieval function
        Branch_Order_fun = @(QSM) QSM.branch.order;

        % Branch orders of initial QSMs
        if ~isempty(QSM_init_cell)
            number_scanner_combinations = length(QSM_init_cell);
            QSM_init_branch_order_cell  = cellfun(Branch_Order_fun, QSM_init_cell, 'UniformOutput', false);
        else
            QSM_init_branch_order_cell  = {};
        end

        % Branch orders of MC QSMs
        if ~isempty(MC_QSM_cell)
            number_scanner_combinations = length(MC_QSM_cell);
            QSM_MC_branch_order_cell    = cell(1, number_scanner_combinations);
    
            for c = 1 : number_scanner_combinations
                QSM_cell                    = MC_QSM_cell{c};
                branch_order_cell           = cellfun(Branch_Order_fun, QSM_cell, 'UniformOutput', false);
                QSM_MC_branch_order_cell{c} = branch_order_cell;
            end
        else
            QSM_MC_branch_order_cell = {};
        end

        % Maximum
        total_branch_order_cell = [QSM_init_branch_order_cell, horzcat(QSM_MC_branch_order_cell{:})];
        branch_order_list       = vertcat(total_branch_order_cell{:});
        max_branch_order        = double(max(branch_order_list));

    %% Branch order parameters %%
        % A structure is created which will contain all the data
        Branch_Order_Parameters     = struct();
        number_parameters           = length(branch_parameters);

        for p = 1 : number_parameters
            % The parameter of interest
            branch_parameter        = branch_parameters{p};
            Branch_Parameter_fun    = @(QSM) QSM.branch.(branch_parameter);

            Branch_Parameter_Data   = struct('Initial', struct(), 'MC', struct());

            for c = 1 : number_scanner_combinations
                % The scanner combination's data
                if ~isempty(QSM_init_cell)
                    QSM_init            = QSM_init_cell{c};
                    init_parameter_list = Branch_Parameter_fun(QSM_init);
                end

                if ~isempty(MC_QSM_cell)
                    MC_QSM_cell_c       = MC_QSM_cell{c};
                    MC_parameter_cell   = cellfun(Branch_Parameter_fun, MC_QSM_cell_c, 'UniformOutput', false);
                end

                % Each branch order's data
                for b = 1 : max_branch_order + 1
                    % Boolean function for this branch order
                    branch_order            = b - 1;
                    Branch_Order_Bool_fun   = @(QSM) QSM.branch.order == branch_order;

                    % Initial QSM
                    if ~isempty(QSM_init_cell)
                        init_branch_order_bool      = Branch_Order_Bool_fun(QSM_init);
                        init_branch_parameter_list  = init_parameter_list(init_branch_order_bool);
    
                        Initial_Parameter_Data = struct('mean', mean(init_branch_parameter_list), 'total', sum(init_branch_parameter_list));
                    else
                        Initial_Parameter_Data = struct();
                    end

                    % Monte Carlo QSMs
                    if ~isempty(MC_QSM_cell)
                        MC_branch_order_bool_cell   = cellfun(Branch_Order_Bool_fun, MC_QSM_cell_c, 'UniformOutput', false);
                        MC_branch_parameter_cell    = Cell_Array_Boolean(MC_parameter_cell, MC_branch_order_bool_cell);
    
                        MC_branch_parameter_mean_list   = cellfun(@mean, MC_branch_parameter_cell);
                        MC_branch_parameter_sum_list    = cellfun(@sum, MC_branch_parameter_cell);
    
                        MC_Parameter_Data = struct('Mean', struct('mean', mean(MC_branch_parameter_mean_list, 'omitnan'), 'std', std(MC_branch_parameter_mean_list, 'omitnan'), 'data', MC_branch_parameter_mean_list), ...
                                                   'Total', struct('mean', mean(MC_branch_parameter_sum_list, 'omitnan'), 'std', std(MC_branch_parameter_sum_list, 'omitnan'), 'data', MC_branch_parameter_sum_list));
                    else
                        MC_Parameter_Data = struct();
                    end

                    % Assigning to the structure
                    Branch_Parameter_Data(c, b) = struct('Initial', Initial_Parameter_Data, 'MC', MC_Parameter_Data);
                end
            end

            % This parameter's data is assigned
            Branch_Order_Parameters.(branch_parameter) = Branch_Parameter_Data;
        end        

    %% Number of branches %%
        % Data
        init_number_branches_matrix     = zeros(number_scanner_combinations, max_branch_order + 1);
        MC_number_branches_mean_matrix  = zeros(number_scanner_combinations, max_branch_order + 1);
        MC_number_branches_std_matrix   = zeros(number_scanner_combinations, max_branch_order + 1);
        MC_number_branches_data_cell    = cell(number_scanner_combinations, max_branch_order + 1);

        for c = 1 : number_scanner_combinations
            % Number of branches per branch order
            for b = 1 : max_branch_order + 1
                branch_order        = b - 1;
                Number_Branches_fun = @(QSM) sum(QSM.branch.order == branch_order);

                % Initial QSM
                if ~isempty(QSM_init_cell)
                    QSM_init                            = QSM_init_cell{c};
                    init_number_branches                = Number_Branches_fun(QSM_init);
                    init_number_branches_matrix(c, b)   = init_number_branches;
                end

                % MC QSMs
                if ~isempty(MC_QSM_cell)
                    MC_QSM_cell_c                           = MC_QSM_cell{c};
                    MC_number_branches_list                 = cellfun(Number_Branches_fun, MC_QSM_cell_c);
                    MC_number_branches_mean_matrix(c, b)    = mean(MC_number_branches_list);
                    MC_number_branches_std_matrix(c, b)     = std(MC_number_branches_list);
                    MC_number_branches_data_cell{c, b}      = MC_number_branches_list;
                end
            end
        end

        % Append to the structure
        Number_Branches = struct('Initial', init_number_branches_matrix, 'MC', struct('mean', MC_number_branches_mean_matrix, 'std', MC_number_branches_std_matrix, 'data', {MC_number_branches_data_cell}));

end