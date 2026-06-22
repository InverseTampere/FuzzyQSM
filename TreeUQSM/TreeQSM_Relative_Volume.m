% The volume of each branch relative to that of all branches is computed in percent and added to the QSM

function [QSM, relative_volume_list] = TreeQSM_Relative_Volume(QSM)

    % Relative volume
    branch_volume_list          = QSM.branch.volume;
    relative_volume_list        = branch_volume_list / sum(branch_volume_list) * 100;
    QSM.branch.RelativeVolume   = relative_volume_list;

end