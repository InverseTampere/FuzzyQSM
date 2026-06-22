% This script determines the cylinder-level, branch-level and tree-level metrics of the smoothed tree model produced by Markku Åkerblom's Blender add-on

function [Smoothed_Model_Metrics, Smoothed_Model_Circle_Geometry, Original_QSM] = Smoothed_Model_Tree_Metrics(QSM_file, smoothed_model_obj_file, Parallel_Pool)

    %% Load the QSM %%
        % Load the QSM
        Original_QSM_File   = load(QSM_file);
        Original_QSM        = Original_QSM_File.QSM;

        % Unchanged by smoothing
        branch_length_list  = Original_QSM.branch.length;       
        branch_order_list   = Original_QSM.branch.order; 

    %% Circles in the object file %%
        % Their geometry
        [Smoothed_Model_Circle_Geometry, Smoothed_Model_Mesh] = Smoothed_Model_Circles(smoothed_model_obj_file);
    
        % Associating them to cylinders
        Smoothed_Model_Circle_Geometry  = Smoothed_Model_Circle_Cylinder_Association(Original_QSM, Smoothed_Model_Circle_Geometry, Parallel_Pool);
        circle_centre_matrix            = Smoothed_Model_Circle_Geometry.centre;
        circle_radius_list              = Smoothed_Model_Circle_Geometry.radius;
        circle_branch_order_list        = Smoothed_Model_Circle_Geometry.BranchOrder;

    %% Segment metrics %%
        Segment_Metrics             = Smoothed_Model_Segment_Metrics(Smoothed_Model_Mesh, Smoothed_Model_Circle_Geometry, Original_QSM, Parallel_Pool);
        segment_volume_list         = Segment_Metrics.volume;
        segment_branch_index_list   = Segment_Metrics.branch_index;
        segment_branch_order_list   = Segment_Metrics.BranchOrder;

    %% Branch-level metrics %%
        % The volume is calculated by summing up the segments within each branch
        number_branches     = length(branch_length_list);
        branch_volume_list  = zeros(1, number_branches);

        for b = 1 : number_branches
            branch_segment_bool         = segment_branch_index_list == b;
            branch_segment_volume_list  = segment_volume_list(branch_segment_bool);
            branch_volume_list(b)       = sum(branch_segment_volume_list);
        end

        % The length and unmodded volume are the original QSM's and regular volume respectively, but kept for consistency
        Branch_Metrics = struct('length', branch_length_list, 'volume', branch_volume_list, 'UnmodVolume', branch_volume_list, 'order', branch_order_list);

    %% Branch-order-level metrics %%
        % Maximum branch order of the smoothed model
        max_branch_order = max(branch_order_list);

        % The lengths and volumes within each branch order are summed
        branch_order_length_list = zeros(1, max_branch_order + 1);
        branch_order_volume_list = zeros(1, max_branch_order + 1);

        for b = 1 : max_branch_order + 1
            % Branches in this order
            branch_order        = b - 1;
            branch_order_bool   = branch_order_list == branch_order;
            
            % Their combined metrics
            branch_order_length         = sum(branch_length_list(branch_order_bool));
            branch_order_length_list(b) = branch_order_length;

            branch_order_volume         = sum(branch_volume_list(branch_order_bool));
            branch_order_volume_list(b) = branch_order_volume;
        end

        % Structure
        Branch_Order_Metrics = struct('length', branch_order_length_list, 'volume', branch_order_volume_list, 'UnmodVolume', branch_order_volume_list, 'max_branch_order', max_branch_order);

    %% Tree metrics %%
        % Total volume
        total_volume = sum(segment_volume_list);

        % Stem volume
        stem_segment_bool   = segment_branch_order_list == 0;
        stem_volume         = sum(segment_volume_list(stem_segment_bool));

        % Branch volume
        branch_volume       = sum(segment_volume_list(~stem_segment_bool));

        % Height
        circle_height_list  = circle_centre_matrix(:, 3);
        tree_height         = max(circle_height_list) - min(circle_height_list);

        % DBH
        stem_circle_bool = circle_branch_order_list == 0;                           % Only consider circles in the stem

        stem_circle_radius_list = circle_radius_list(stem_circle_bool);
        stem_circle_height_list = circle_height_list(stem_circle_bool);
        [~, circle_order]       = sort(stem_circle_height_list);

        BH              = 1.3;
        circle_above_BH = circle_order(find(stem_circle_height_list > BH, 1));      % The first circle above BH
        circle_below_BH = circle_above_BH - 1;                                      % The last circle below it

        height_above = stem_circle_height_list(circle_above_BH);
        radius_above = stem_circle_radius_list(circle_above_BH);
        
        height_below = stem_circle_height_list(circle_below_BH);
        radius_below = stem_circle_radius_list(circle_below_BH);

        weight_above = 1 - (height_above - BH) / (height_above - height_below);     % Weighted inversely to the difference to breast-height
        weight_below = 1 - (BH - height_below) / (height_above - height_below);

        DBH = 2*(weight_above*radius_above + weight_below*radius_below);

        % Structure
        Tree_Metrics = struct('TotalVolume', total_volume, 'TrunkVolume', stem_volume, 'BranchVolume', branch_volume, 'TotalUnmodVolume', total_volume, 'TrunkUnmodVolume', stem_volume, 'BranchUnmodVolume', branch_volume, 'DBHqsm', DBH, 'UnmodDBHqsm', DBH, 'TreeHeight', tree_height);

    %% Structure containing all the metrics %%
        Smoothed_Model_Metrics = struct('treedata', Tree_Metrics, 'branch', Branch_Metrics, 'BranchOrder', Branch_Order_Metrics, 'segment', Segment_Metrics);

end