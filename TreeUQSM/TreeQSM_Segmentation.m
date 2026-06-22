% This file is a modification of TreeQSM where it only creates cover sets and segments
% This then separates this stochastic element from the deterministic cylinder fitting
% Note that only one input should be given, i.e. the optimal after using select_optimum

function [Cover_Sets, Segments] = TreeQSM_Segmentation(P, Opt_Inputs)

%% Ensure the point cloud is the proper form %%
    % only 3-dimensional data
    if size(P,2) > 3
        P = P(:,1:3);
    end
    % Only double precision data
    if ~isa(P,'double')
        P = double(P);
    end
  
  %% Generate cover sets %%
      cover1 = cover_sets(P, Opt_Inputs);
  
  %% Determine tree sets and update neighbors %%
    [cover1, Base, Forb] = tree_sets(P, cover1, Opt_Inputs);
  
  %% Determine initial segments
      segment1 = segments(cover1, Base, Forb);
  
  %% Correct segments %%
      % Don't remove small segments and add the modified base to the segment
      segment1 = correct_segments(P, cover1, segment1, Opt_Inputs, 0, 1, 1);

  %% Generate new cover sets %%
      % Determine relative size of new cover sets and use only tree points
      RS = relative_size(P, cover1, segment1);
      
      % Generate new cover
      Cover_Sets = cover_sets(P, Opt_Inputs, RS);
      
  %% Determine tree sets and update neighbors %%
      [Cover_Sets, Base, Forb] = tree_sets(P, Cover_Sets, Opt_Inputs, segment1);
      
  %% Determine segments %%
      Segments = segments(Cover_Sets, Base, Forb);
      
  %% Correct segments
      % Remove small segments and the extended bases.
      Segments = correct_segments(P, Cover_Sets, Segments, Opt_Inputs, 1, 1, 0);
end
