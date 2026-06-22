% This script retrieves the point clouds and other properties of Helios generated .xyz files as well as the scanner positions in the associated .txt file
% The files in the given folder are expected to be of the standard format, i.e. leg000_points.xyz and leg000_trajectory.txt

% The output fields of Point_Cloud_Data:
%       point_cloud_cell            n x 3
%       number_points_list (n)      1 x 1
%       intensity_cell              n x 1
%       echo_width_cell             n x 1
%       return_number_cell          n x 1
%       number_of_returns_cell      n x 1
%       full_wave_index_cell        n x 1
%       hit_object_ID_cell          n x 1
%       class_cell                  n x 1
%       GPS_time_cell               n x 1

function Point_Cloud_Data = Helios_Data_Reader(folder_path, coarsening_factor, Plot)

    %% Retrieve the data %%
        % Point cloud files
        point_cloud_files   = dir(sprintf('%s/leg*_points.xyz', folder_path));
        number_scanners     = length(point_cloud_files);

        % Cell arrays for retrieved data
        [point_cloud_cell, intensity_cell, echo_width_cell, return_number_cell, number_returns_cell, full_wave_index_cell, hit_object_ID_cell, class_cell, GPS_time_cell] = deal(cell(1, number_scanners));

        for s = 1 : number_scanners
            % Lines of data in the file
            point_cloud_file_name   = point_cloud_files(s).name;
            point_cloud_file_path   = sprintf('%s/%s', folder_path, point_cloud_file_name);
            point_cloud_data_matrix = readmatrix(point_cloud_file_path, 'FileType', 'text');

            % Coarsened
            if coarsening_factor > 1
                point_cloud_data_matrix = point_cloud_data_matrix(1 : coarsening_factor : end, :);
            end
    
            % Data types
            point_cloud_cell{s}     = point_cloud_data_matrix(:, 1 : 3);
            intensity_cell{s}       = point_cloud_data_matrix(:, 4);
            echo_width_cell{s}      = point_cloud_data_matrix(:, 5);
            return_number_cell{s}   = point_cloud_data_matrix(:, 6);
            number_returns_cell{s}  = point_cloud_data_matrix(:, 7);
            full_wave_index_cell{s} = point_cloud_data_matrix(:, 8);
            hit_object_ID_cell{s}   = point_cloud_data_matrix(:, 9);
            class_cell{s}           = point_cloud_data_matrix(:, 10);
            GPS_time_cell{s}        = point_cloud_data_matrix(:, 11);
        end

        % Structure
        number_points_list  = cellfun(@length, point_cloud_cell);
        Point_Cloud_Data    = struct('point_cloud_cell', {point_cloud_cell}, 'number_points_list', number_points_list, 'intensity_cell', {intensity_cell}, 'echo_width_cell', {echo_width_cell}, 'return_number_cell', {return_number_cell}, 'number_returns_cell', {number_returns_cell}, 'full_wave_index_cell', {full_wave_index_cell}, 'hit_object_ID_cell', {hit_object_ID_cell}, 'class_cell', {class_cell}, 'GPS_time_cell', {GPS_time_cell});

    %% Plots %%
        if Plot == true
            % Colour map
            point_cloud_cmap = cbrewer('qual', 'Set1', max(number_scanners, 3));

            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])  
            
            hold on
            grid on

            % Each scanner's location and point cloud
            for s = 1 : number_scanners
                point_cloud_matrix  = point_cloud_cell{s};
                scanner_colour      = point_cloud_cmap(s, :);
    
                % Point cloud
                scatter3(point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 5, 'MarkerFaceColor', scanner_colour, 'MarkerEdgeColor', 'none', 'DisplayName', sprintf('Scanner %i', s));
            end
            
            % Axes
            xlabel('x [m]');
            ylabel('y [m]');
            zlabel('z [m]');
    
            axis equal
            view(45, 45);

            % Legend
            legend('show', 'location', 'eastoutside');
    
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
    
            hold off
    
            % Pause message
            disp('The point cloud has been retrieved. The script will continue and figure close upon a key-press.');
            pause();
            
            close(1);
        end
end