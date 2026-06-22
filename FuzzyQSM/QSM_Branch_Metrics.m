% This script determines branch and branch-order level length and volume using the reference QSMs topology

function [Branch_Metrics, Branch_Order_Metrics] = QSM_Branch_Metrics(Cylinder_Geometry, QSM)

    %% Structure inputs %%
        % QSM
        cylinder_branch_index_list  = QSM.cylinder.branch;
        branch_order_list           = QSM.branch.order;

        % Cylinder geometry
        cylinder_length_list        = Cylinder_Geometry.length;
        cylinder_radius_list        = Cylinder_Geometry.radius;

        if isfield(Cylinder_Geometry, 'UnmodRadius')                        % It may not be present
            cylinder_UnmodRadius_list   = Cylinder_Geometry.UnmodRadius;
        else
            cylinder_UnmodRadius_list   = Cylinder_Geometry.radius;
        end

    %% Branch metrics %%
        % Length and volume
        cylinder_volume_list        = pi*cylinder_radius_list.^2 .* cylinder_length_list;
        cylinder_UnmodVolume_list   = pi*cylinder_UnmodRadius_list.^2 .* cylinder_length_list;

        number_branches         = length(branch_order_list);
        branch_length_list      = zeros(1, number_branches);
        branch_volume_list      = zeros(1, number_branches);
        branch_UnmodVolume_list = zeros(1, number_branches);

        for b = 1 : number_branches
            % This branch' cylinders
            branch_cylinder_bool        = cylinder_branch_index_list == b;

            branch_cylinder_length_list = cylinder_length_list(branch_cylinder_bool);
            branch_length_list(b)       = sum(branch_cylinder_length_list);

            branch_cylinder_volume_list = cylinder_volume_list(branch_cylinder_bool);
            branch_volume_list(b)       = sum(branch_cylinder_volume_list);

            branch_cylinder_UnModVolume_list    = cylinder_UnmodVolume_list(branch_cylinder_bool);
            branch_UnmodVolume_list(b)          = sum(branch_cylinder_UnModVolume_list);
        end

        % A structure is created
        Branch_Metrics = struct('length', branch_length_list, 'volume', branch_volume_list, 'UnmodVolume', branch_UnmodVolume_list, 'order', branch_order_list);

    %% Branch-order metrics %%
        % Summed values for each branch order
        max_branch_order = max(branch_order_list);

        branch_order_length_list        = zeros(1, max_branch_order + 1);
        branch_order_volume_list        = zeros(1, max_branch_order + 1);
        branch_order_UnmodVolume_list   = zeros(1, max_branch_order + 1);

        for b = 1 : max_branch_order + 1
            % Branches for this branch order
            branch_order        = b - 1;
            branch_order_bool   = branch_order_list == branch_order;

            % Their metrics
            branch_order_length         = branch_length_list(branch_order_bool);
            branch_order_length_list(b) = sum(branch_order_length);
            
            branch_order_volume         = branch_volume_list(branch_order_bool);
            branch_order_volume_list(b) = sum(branch_order_volume);

            branch_order_UnmodVolume            = branch_UnmodVolume_list(branch_order_bool);
            branch_order_UnmodVolume_list(b)    = sum(branch_order_UnmodVolume);
        end

        % Structure
        Branch_Order_Metrics = struct('length', branch_order_length_list, 'volume', branch_order_volume_list, 'UnmodVolume', branch_order_UnmodVolume_list, 'max_branch_order', max_branch_order);

end