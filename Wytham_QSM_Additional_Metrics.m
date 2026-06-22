% This script goes through the original QSMs and those fitted by TreeQSM and FuzzyQSM to add new metrics based on existing information
% Note that the script goes through all QSMs ending in QSM.mat located in the given folders

clear variables
close all
clc

%% Inputs %%
    % Data folders
    Original_QSMs_folder    = 'Data/Original_QSMs';
    TreeQSMs_folder         = 'Data/TreeQSMs';
    FuzzyQSMs_folder        = 'Data/FuzzyQSMs';

%% Appending to the QSMs %%
    % Cell array containing all folders
    QSM_type_folders = {Original_QSMs_folder, TreeQSMs_folder, FuzzyQSMs_folder};
    number_QSM_types = length(QSM_type_folders);

    for t = 1 : number_QSM_types
        % All QSMs located in this folder
        QSM_type_folder = QSM_type_folders{t};

        QSM_Files   = dir(sprintf('%s/**/*QSM.mat', QSM_type_folder));
        number_QSMs = length(QSM_Files);

        for q = 1 : number_QSMs
            % Load the QSM
            QSM_File    = QSM_Files(q);
            QSM_folder  = QSM_File.folder;
            QSM_name    = QSM_File.name;

            QSM_file_path   = sprintf("%s%s%s", QSM_folder, '\', QSM_name);
            QSM_Data        = load(QSM_file_path);

            if contains(QSM_type_folder, 'Original')
                QSM = QSM_Data.QSM;
            elseif contains(QSM_type_folder, 'Discrete')
                QSM = QSM_Data.QSM_Discrete;
            elseif contains(QSM_type_folder, 'Fuzzy')
                QSM = QSM_Data.QSM_Fuzzy;
            end

            % Apply the functions
            QSM = Total_Branch_Length(QSM);

            % Save it
            if contains(QSM_type_folder, 'Original')
                QSM_Data.QSM = QSM;
            elseif contains(QSM_type_folder, 'Discrete')
                QSM_Data.QSM_Discrete = QSM;
            elseif contains(QSM_type_folder, 'Fuzzy')
                QSM_Data.QSM_Fuzzy = QSM;
            end

            save(QSM_file_path, '-struct', 'QSM_Data');
        end
    end

%% Applied functions %%
    function QSM = Total_Branch_Length(QSM)
        % Total branch length
        branch_length_list  = QSM.branch.length;
        total_branch_length = sum(branch_length_list);

        % Added to the tree metrics
        QSM.treedata.BranchLength = total_branch_length;
    end