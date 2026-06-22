% This script compensates for the range bias within the data

function point_cloud_cell = Range_Bias_Compensation_3D(point_cloud_cell, range_bias, Scanner_loc_cell)

    %% Range bias compensation for each scanner %%
        number_scanners = length(Scanner_loc_cell);

        if range_bias ~= 0
            for s = 1 : number_scanners
                % The point cloud is moved closer or farther away from the scanner
                scanner_loc     = Scanner_loc_cell{s};
                points_matrix   = point_cloud_cell{s};

                % The normalised vectors from the scanner to the points
                vector_matrix   = points_matrix - scanner_loc;

                vector_matrix_n = vector_matrix ./ sqrt(sum(vector_matrix.^2, 2));

                % The range bias is compensated for
                bias_comp_matrix    = range_bias * vector_matrix_n;
                points_matrix_comp  = points_matrix - bias_comp_matrix;

                point_cloud_cell{s} = points_matrix_comp;
            end
        end
end