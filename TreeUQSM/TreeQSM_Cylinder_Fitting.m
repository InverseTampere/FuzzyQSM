% This file follows after TreeQSM Segmentation. Cylinder fitting is deterministic, whilst segmentation is stochastic
% Note that only one input should be given, i.e. the optimal after using select_optimum

function QSM = TreeQSM_Cylinder_Fitting(P, Cover_Sets, Segments, Opt_Inputs)
      
      %% Define cylinders %%
          cylinder = cylinders(P, Cover_Sets, Segments, Opt_Inputs);
      
      if ~isempty(cylinder.radius)
        %% Determine the branches
            branch = branches(cylinder);
        
        %% Compute (and display) model attributes %%
            T = Segments.segments{1};
            T = vertcat(T{:});
            T = vertcat(Cover_Sets.ball{T});
            trunk = P(T,:); % point cloud of the trunk
            % Compute attributes and distibutions from the cylinder model
            % and possibly some from a triangulation
            [treedata,triangulation] = tree_data(cylinder,branch,trunk,Opt_Inputs);
        
        %% Compute point model distances %%
            if Opt_Inputs.Dist
              pmdis = point_model_distance(P,cylinder);
            end
        
        %% Reconstruct the output "QSM" %%
            clear QSM
            QSM = struct('cylinder',{},'branch',{},'treedata',{},'rundata',{},...
              'pmdistance',{},'triangulation',{});
            QSM(1).cylinder = cylinder;
            QSM(1).branch = branch;
            QSM(1).treedata = treedata;
            QSM(1).rundata.Opt_Inputs = Opt_Inputs;
            QSM(1).rundata.version = '2.4.1';
        
            if Opt_Inputs.Dist
              QSM(1).pmdistance = pmdis;
            end
            if Opt_Inputs.Tria
              QSM(1).triangulation = triangulation;
            end
        
        %% Save the output into results-folder %%
            % matlab-format (.mat)
            if Opt_Inputs.savemat
              str = [Opt_Inputs.name,'_t',num2str(Opt_Inputs.tree),'_m',...
                num2str(Opt_Inputs.model)];
              save(['results/QSM_',str],'QSM')
            end

        % text-format (.txt)
        if Opt_Inputs.savetxt
            str = [Opt_Inputs.name,'_t',num2str(Opt_Inputs.tree),'_m',...
            num2str(Opt_Inputs.model)];
            save_model_text(QSM,str)
        end

        %% Plot models and segmentations %%
            if Opt_Inputs.plot >= 1
                if Opt_Inputs.Tria
                    plot_models_segmentations(P,Cover_Sets,Segments,cylinder,trunk,...
                    triangulation)
                else
                    plot_models_segmentations(P,Cover_Sets,Segments,cylinder)
                end
            end
      end
end