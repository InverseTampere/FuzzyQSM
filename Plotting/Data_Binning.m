% This function determines the average of the given data for the given bins in x

function [data_avg_list, number_points_bin_list, bin_bool_cell] = Data_Binning(data_list, x_list, x_bin_edges_list, number_bins)
    % The average value within each bin as well as the (number of) members within each bin are saved
    data_avg_list           = zeros(number_bins, 1);
    number_points_bin_list  = zeros(number_bins, 1);
    bin_bool_cell           = cell(1, number_bins);
    
    for b = 1 : number_bins
        % The bin edges
        x_lb = x_bin_edges_list(b);
        x_ub = x_bin_edges_list(b + 1);

        % Entries within the bin
        bin_bool            = x_list > x_lb & x_list < x_ub;
        bin_bool_cell{b}    = bin_bool;
        
        % The average
        bin_data            = data_list(bin_bool);
        data_avg_list(b)    = mean(bin_data);

        % The number of points within the bin
        number_points               = sum(bin_bool);
        number_points_bin_list(b)   = number_points;
    end
end