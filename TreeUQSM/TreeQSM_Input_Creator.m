% Creates input parameter structure array needed to run TreeQSM function

function TreeQSM_Inputs = TreeQSM_Input_Creator(Fitting_Parameters, Data_Parameters, PG_Inputs_Defined, point_cloud_matrix, tree_identifier, TreeQSM_Print, TreeQSM_Plot)

    %% Inputs %%
        % Fitting parameters
        patch_diameter_1            = Fitting_Parameters.patch_diameter_1;
        patch_diameter_2_min        = Fitting_Parameters.patch_diameter_2_min;
        patch_diameter_2_max        = Fitting_Parameters.patch_diameter_2_max;

        ball_radius_1_increase      = Fitting_Parameters.ball_radius_1_increase;
        min_number_ball_points_1    = Fitting_Parameters.min_number_ball_points_1;
        ball_radius_2_increase      = Fitting_Parameters.ball_radius_2_increase;
        min_number_ball_points_2    = Fitting_Parameters.min_number_ball_points_2;

        min_cylinder_radius         = Fitting_Parameters.min_cylinder_radius;
        Parent_correction           = Fitting_Parameters.Parent_correction;
        Taper_correction            = Fitting_Parameters.Taper_correction;

        Growth_volume_correction    = Fitting_Parameters.Growth_volume_correction;
        growth_volume_factor        = Fitting_Parameters.growth_volume_factor;

        % Data parameters
        Only_tree                   = Data_Parameters.Only_tree;

    %% Creating the inputs structure %%
        % Patch generation parameters
        if PG_Inputs_Defined == true
            % Number of desired inputs
            [nPD1, nPD2Min, nPD2Max] = deal(length(patch_diameter_1), length(patch_diameter_2_min), length(patch_diameter_2_max));

            % Automatically defined inputs for patch generation
            Defined_Inputs = define_input(point_cloud_matrix, nPD1, nPD2Min, nPD2Max);

            [patch_diameter_1, patch_diameter_2_min, patch_diameter_2_max]                      = deal(Defined_Inputs.PatchDiam1, Defined_Inputs.PatchDiam2Min, Defined_Inputs.PatchDiam2Max);
            [ball_radius_1, min_number_ball_points_1, ball_radius_2, min_number_ball_points_2]  = deal(Defined_Inputs.BallRad1, Defined_Inputs.nmin1, Defined_Inputs.BallRad2, Defined_Inputs.nmin2);
            [ball_radius_1_increase, ball_radius_2_increase]                                    = deal(ball_radius_1 - patch_diameter_1, ball_radius_2 - patch_diameter_2_max);
        end

        TreeQSM_Inputs.PatchDiam1       = patch_diameter_1; 
        TreeQSM_Inputs.PatchDiam2Min    = patch_diameter_2_min; 
        TreeQSM_Inputs.PatchDiam2Max    = patch_diameter_2_max; 
    
        TreeQSM_Inputs.BallRad1     = patch_diameter_1 + ball_radius_1_increase; 
        TreeQSM_Inputs.nmin1        = min_number_ball_points_1; 
        TreeQSM_Inputs.BallRad2     = patch_diameter_2_max + ball_radius_2_increase; 
        TreeQSM_Inputs.nmin2        = min_number_ball_points_2; 

        % Radius correction
        TreeQSM_Inputs.MinCylRad    = min_cylinder_radius; 
        TreeQSM_Inputs.ParentCor    = Parent_correction; 
        TreeQSM_Inputs.TaperCor     = Taper_correction; 

        % Growth volume correction
        TreeQSM_Inputs.GrowthVolCor = Growth_volume_correction;
        TreeQSM_Inputs.GrowthVolFac = growth_volume_factor;

        % Other inputs
        TreeQSM_Inputs.OnlyTree     = double(Only_tree); 
        TreeQSM_Inputs.name         = tree_identifier; 

        TreeQSM_Inputs.plot         = 2*double(TreeQSM_Plot);
        TreeQSM_Inputs.disp         = 2*double(TreeQSM_Print);

        TreeQSM_Inputs.Tria         = 0;                        % Not currently required
        TreeQSM_Inputs.Dist         = 1;                        % Required
        TreeQSM_Inputs.tree         = 1;                        % Tree indexing
        TreeQSM_Inputs.model        = 1;                        % Model indexing
        TreeQSM_Inputs.savemat      = 1;                        % Potentially used
        TreeQSM_Inputs.savetxt      = 0;                        % Not currently used

end
