% The total volume of the cylinders within the given height bins is determined
% The centre of the cylinders is used and the zeroth order cylinders (tree stem) are split from the branches

% Note that the heights are adjusted s.t. the lowest cylinder has a centre at 0

function [branch_total_volume_list, branch_number_cylinders_list, stem_total_volume_list, stem_number_cylinders_list] = QSM_Height_Volume(QSM, bin_edge_heights)

    %% Inputs %%
        % QSM
        cylinder_start_matrix   = QSM.cylinder.start;
        cylinder_axis_matrix    = QSM.cylinder.axis;
        cylinder_length_list    = QSM.cylinder.length;
        cylinder_radius_list    = QSM.cylinder.radius;
        branch_order_list       = QSM.cylinder.BranchOrder;

    %% Binning %%
        % Volume of each cylinder  
        cylinder_volume_list    = pi*cylinder_radius_list.^2 .* cylinder_length_list;

        % Height of each cylinder
        cylinder_centre_matrix  = cylinder_start_matrix + 1/2*cylinder_length_list .* cylinder_axis_matrix;
        cylinder_height_list    = cylinder_centre_matrix(:, 3);
        cylinder_height_list    = cylinder_height_list - min(cylinder_height_list);

        % Tree stem cylinders have a branch order of zero
        tree_stem_bool = branch_order_list == 0;

        % Total volume within each bin
        number_bins                     = length(bin_edge_heights) - 1;
        branch_total_volume_list        = zeros(1, number_bins);
        branch_number_cylinders_list    = zeros(1, number_bins);
        stem_total_volume_list          = zeros(1, number_bins);
        stem_number_cylinders_list      = zeros(1, number_bins);

        for i = 1 : number_bins
            % Height range of the bin
            height_bottom   = bin_edge_heights(i);
            height_top      = bin_edge_heights(i + 1);

            % Cylinders within this bin
            height_bin_bool = cylinder_height_list > height_bottom & cylinder_height_list <= height_top;
            branch_bin_bool = height_bin_bool & ~tree_stem_bool;
            stem_bin_bool   = height_bin_bool & tree_stem_bool;

            % Volume
            branch_total_volume_list(i) = sum(cylinder_volume_list(branch_bin_bool));
            stem_total_volume_list(i)   = sum(cylinder_volume_list(stem_bin_bool));

            % Number of cylinders
            branch_number_cylinders_list(i) = sum(branch_bin_bool);
            stem_number_cylinders_list(i)   = sum(stem_bin_bool);
        end
end