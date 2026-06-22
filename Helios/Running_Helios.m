% This script runs Helios for the given survey and scene files by editing them to scan the given tree with the given settings (beam divergence and scanner locations)

% Inputs:
%   Helios_folder                                   Folder in which Helios is installed, i.e. "C:\Users\svvive\AppData\Local\helios"
%   scene_file                                      Scene file path within the Helios folder, i.e. "data\scenes\Wytham_Woods\Wytham_Scene.xml"
%   survey_file                                     Survey file path within the Helios folder, i.e. "data\surveys\Wytham_Woods\Wytham_Survey.xml"
%   scanner_file                                    TLS scanners file path within the Helios folder, i.e. "Lib\site-packages\pyhelios\data\scanners_tls.xml"
%   tree_file                                       Tree .obj file path, i.e. "C:\Users\svvive\OneDrive - TUNI.fi\Remote_Sensing_Forests\FuzzyQSM\Data\Smoothed_Models\ACERPS\174a.mtl" 
%   output_folder                                   Folder into which the output is moved and renamed, i.e. "C:\Users\svvive\OneDrive - TUNI.fi\Remote_Sensing_Forests\FuzzyQSM\Data\Helios_Scans"
%   Scanning_Parameters:
%       beam_divergence                             Beam divergence full-angle in radians
%       Scanner_loc_cell            1 x 3           Note that the final scanner locations include the tripod
%       number_scanners             1 x 1

% Outputs:
%   Point_Cloud_Data:
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
%   Scanning_Parameters:
%       Scanner_loc_cell            1 x 3
%       number_scanners             1 x 1

function Point_Cloud_Data = Running_Helios(Helios_folder, scene_file, survey_file, scanner_file, scanner_ID, tree_file, Helios_batch_file, output_folder, Scanning_Parameters, Plot)
    
    %% Structure inputs %%
        % Scanning parameters
        beam_divergence     = Scanning_Parameters.beam_divergence;
        Scanner_loc_cell    = Scanning_Parameters.Scanner_loc_cell;
        number_scanners     = Scanning_Parameters.number_scanners;

    %% Edit the scanner file %%
        % Scanner file lines
        scanner_file_path   = sprintf('%s%s%s', Helios_folder, '\', scanner_file);
        scanner_file_lines  = readlines(scanner_file_path);
        Scanner_Finder      = @(line) contains(line, scanner_ID);                                                                   % Find where the desired scanner is specified
        scanner_start_index = find(cellfun(Scanner_Finder, scanner_file_lines));                
        
        Beam_Divergence_Finder      = @(line) contains(line, 'beamDivergence_rad');                                                 % Find lines specifying the beam divergence
        beam_divergence_index_list  = find(cellfun(Beam_Divergence_Finder, scanner_file_lines));
        beam_divergence_index       = beam_divergence_index_list(find(beam_divergence_index_list > scanner_start_index, 1));        % This scanner's beam divergence line

        % Update the beam divergence and scanner file
        beam_divergence_line                        = sprintf("            beamDivergence_rad         = ""%f""", beam_divergence);
        scanner_file_lines(beam_divergence_index)   = beam_divergence_line;
        
        writelines(scanner_file_lines, scanner_file_path);

    %% Edit the scene file %%
        % Scene file lines
        scene_file_path     = sprintf('%s%s%s', Helios_folder, '\', scene_file);
        scene_file_lines    = readlines(scene_file_path);
        Scene_Start_Finder  = @(line) contains(line, '<scene');                     % Find the beginning of the scene       
        scene_start_index   = find(cellfun(Scene_Start_Finder, scene_file_lines));
        Scene_End_Finder    = @(line) contains(line, '</scene>');                   % Find the end of the scene       
        scene_end_index     = find(cellfun(Scene_End_Finder, scene_file_lines));

        %--% Scanned tree %--%
        % Lines to define the tree object part
        tree_object_ID    = 1;
        tree_object_lines = [sprintf("\t <part id=""%i"">", tree_object_ID);
                             sprintf("\t \t <filter type=""objloader"">");
                             sprintf("\t \t \t <param type=""string"" key=""filepath"" value=""%s""/>", tree_file);
				             sprintf("\t \t \t <param type=""string"" key=""up"" value=""y"" />");                      % y-axis is defined as up due to how Blender defines objects
                             sprintf("\t \t </filter>");
			                 sprintf("\t \t <filter type=""translate"">");
                             sprintf("\t \t \t <param type=""vec3"" key=""offset"" value=""0;0;0"" />");
                             sprintf("\t \t </filter>");
                             sprintf("\t </part>")];

        %--% Complete scene file %--%
        % .xml file
        updated_scene_lines = [scene_file_lines(1 : scene_start_index); tree_object_lines; scene_file_lines(scene_end_index : end)];
        writelines(updated_scene_lines, scene_file_path);

        % .scene file is deleted if it exists
        scene_file_name = strrep(scene_file_path, '.xml', '.scene');

        if exist(scene_file_name, 'file')
            delete(scene_file_name);
        end

    %% Edit the survey file %%
        % Survey file lines
        survey_file_path    = sprintf('%s%s%s', Helios_folder, '\', survey_file);
        survey_file_lines   = readlines(survey_file_path);

        % Lines beginning and ending the survey
        Survey_Start_Finder = @(line) contains(line, '<survey');
        survey_start_index  = find(cellfun(Survey_Start_Finder, survey_file_lines), 1);
        Survey_End_Finder   = @(line) contains(line, '</survey');
        survey_end_index    = find(cellfun(Survey_End_Finder, survey_file_lines), 1);

        % Specifying the correct scanner
        survey_start_line   = survey_file_lines{survey_start_index};
        survey_start_line   = strrep(survey_start_line, '<', '');                       % Remove <> at beginning and end
        survey_start_line   = strrep(survey_start_line, '>', '');                       
        survey_line_parts   = strsplit(survey_start_line, ' ');
        Scanner_Finder      = @(part) contains(part, 'scanner');
        scanner_bool        = arrayfun(Scanner_Finder, survey_line_parts);
        scanner_string      = sprintf('scanner="%s#%s"', scanner_file, scanner_ID);
        
        survey_line_parts{scanner_bool} = scanner_string;
        survey_start_line               = strjoin(survey_line_parts, ' ');
        survey_start_line               = sprintf("%s>", survey_start_line);           % Add <> at beginning and end again
        survey_start_line               = strrep(survey_start_line, 'survey', '<survey');

        % Optional line containing the full waveform settings
        FWF_Line_Finder     = @(line) contains(line, 'FWFSettings');
        FWF_line_bool       = cellfun(FWF_Line_Finder, survey_file_lines);

        % The scanner location code
        scanner_leg_cell = cell(1, number_scanners);

        for s = 1 : number_scanners
            scanner_loc = Scanner_loc_cell{s};

            scanner_leg_lines   = [sprintf("\t \t <leg>");
                                   sprintf("\t \t \t <platformSettings x=""%f"" y=""%f"" z=""%f"" />", scanner_loc(1), scanner_loc(2), scanner_loc(3));
                                   sprintf("\t \t \t <scannerSettings template=""Scanner_Profile"" />");
                                   sprintf("\t \t </leg>")];
            scanner_leg_cell{s} = scanner_leg_lines;
        end

        scanner_legs_lines = vertcat(scanner_leg_cell{:});

        % Updated survey file
        updated_survey_lines = [survey_file_lines(1 : survey_start_index - 1); survey_start_line; survey_file_lines(FWF_line_bool); scanner_legs_lines; survey_file_lines(survey_end_index : end)];
        writelines(updated_survey_lines, survey_file_path);

    %% Edit the batch script %%
        % Lines of code in the batch file
        Helios_batch_file_path  = sprintf('%s%s%s', Helios_folder, '\', Helios_batch_file);
        batch_file_lines        = readlines(Helios_batch_file_path);
    
        % Find the executive line that runs Helios
        Exec_Line_Finder    = @(line) contains(line, '"%windir%\system32\cmd.exe" /C ""{{ ROOT_PREFIX }}\Scripts\activate.bat"');       
        exec_line_bool      = arrayfun(Exec_Line_Finder, batch_file_lines);
        exec_line           = batch_file_lines{exec_line_bool};
        exec_line           = exec_line(2 : end - 1);                       % Removes " at start and end
    
        % Split it to find the Helios command within it
        exec_line_parts         = strsplit(exec_line, '&');                                                             % Ampersand is used to parse commands
        Helios_Command_Finder   = @(part) strcmp(part(2:7), 'helios') & strcmp(part(end - 3 : end), '.xml');            % To ensure that the right command is edited it starts with ' helios' and ends with the survey file extension
        Helios_command_bool     = cellfun(Helios_Command_Finder, exec_line_parts);
    
        % Edit it for the new survey file
        Helios_command                          = sprintf(' helios %s', survey_file_path);      % New Helios command
        exec_line_parts{Helios_command_bool}    = Helios_command;                               % Substitutes the Helios command
        exec_line                               = strjoin(exec_line_parts, '&');                % Joins the parts together
        exec_line                               = sprintf('"%s"', exec_line);                   % Adds the " at start and end
        
        % Write the new batch script
        batch_file_lines{exec_line_bool} = exec_line;
        writelines(batch_file_lines, Helios_batch_file_path);

    %% Run Helios %%
        % Run the batch script
        [status, cmdout] = system(sprintf('cmd /c ""%s""', Helios_batch_file_path), '-echo');  

        if status ~= 0
            error("Running the Helios file failed with status %d.\nOutput:\n%s", status, cmdout);
        end

    %% Read the data %%
        % Folder in which the outputs are placed
        Survey_ID_Finder        = @(part) contains(part, 'name');
        survey_ID_bool          = cellfun(Survey_ID_Finder, survey_line_parts);
        survey_ID_part          = survey_line_parts{survey_ID_bool};
        survey_ID               = strrep(survey_ID_part, '=', '');
        survey_ID               = strrep(survey_ID, 'name', '');
        survey_ID               = strrep(survey_ID, '"', '');

        survey_ID_output_folder         = sprintf('%s%soutput%s%s', Helios_folder, '\', '\', survey_ID);                % Contains all results for this survey ID
        survey_ID_output_folder_files   = dir(survey_ID_output_folder);
        survey_output_folder            = survey_ID_output_folder_files(end).name;                                      % It should always be the last entry
        Helios_output_folder            = sprintf("%s%s%s", survey_ID_output_folder, '\', survey_output_folder);

        % The data is read but not coarsened
        coarsening_factor   = 1;
        Point_Cloud_Data    = Helios_Data_Reader(Helios_output_folder, coarsening_factor, Plot);

        % The folder is moved and renamed
        tree_file_name_parts    = strsplit(tree_file, '\');
        tree_file_obj           = tree_file_name_parts{end};
        tree_ID                 = strrep(tree_file_obj, '.obj', '');
        output_folder_name      = sprintf("%s%s%s", output_folder, '\', tree_ID);

        movefile(Helios_output_folder, output_folder_name);
        
end