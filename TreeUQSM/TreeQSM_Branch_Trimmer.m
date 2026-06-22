% TreeQSM includes branches in the structure which do not actually exist and are thus removed
% They are recognised by the area

function QSM = TreeQSM_Branch_Trimmer(QSM)

    %% Structure inputs %%
        % Branches
        Branches                    = QSM.branch;
        branch_area_list            = Branches.area;

        % Cylinders
        cylinder_branch_index_list  = QSM.cylinder.branch;
    
    %% Remove empty branches %%
        % Only branches with an area are kept
        existant_bool   = branch_area_list > 0;
        Branches        = Structure_Boolean(Branches, existant_bool);
        QSM.branch      = Branches;

        branch_parent_list = Branches.parent;

    %% Change branch indices %%
        % Unique branches
        unique_branches = unique(cylinder_branch_index_list);           % Note that this also automatically sorts
        number_branches = length(unique_branches);

        % Assign new branch indices
        number_cylinders            = length(cylinder_branch_index_list);
        new_cyl_branch_index_list   = zeros(number_cylinders, 1);
        new_branch_parent_list      = zeros(number_branches, 1);

        for b = 1 : number_branches
            % This branch's old index
            old_branch_index = unique_branches(b);

            % Change the cylinder branch indices
            cyl_branch_bool                             = cylinder_branch_index_list == old_branch_index;
            new_cyl_branch_index_list(cyl_branch_bool)  = b;

            % Change the branch parent indices
            branch_parent_bool                          = branch_parent_list == old_branch_index;
            new_branch_parent_list(branch_parent_bool)  = b;
        end

        % Updating the structure
        QSM.cylinder.branch = new_cyl_branch_index_list;
        QSM.branch.parent   = new_branch_parent_list;

end