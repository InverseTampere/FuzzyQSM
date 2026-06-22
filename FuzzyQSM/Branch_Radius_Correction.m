% Taper and parent correction are applied to the branch cylinders if desired akin to TreeQSM
% Note that most of it is taken from adjustments in cylinders.m of TreeQSM v2.4.1

function Cylinders = Branch_Radius_Correction(Cylinders, TreeQSM_Inputs, cyl_point_cloud_cell)

    %% Inputs %%
        % Cylinder data
        cylinder_branch_list    = Cylinders.branch;
        cylinder_parent_list    = Cylinders.parent;
        cylinder_position_list  = Cylinders.PositionInBranch;
        cylinder_BO_list        = Cylinders.BranchOrder;
        cylinder_radius_list    = Cylinders.radius;
        cylinder_length_list    = Cylinders.length;
        cylinder_centre_matrix  = Cylinders.centre;
        cylinder_axis_matrix    = Cylinders.axis;

    %% Inputs %%
        Diagnostics     = false;             % [true, false] Shows the result for a branch

    %% Surface coverage %%
        % Same inputs as TreeQSM's least_squares_cylinder
        res = 0.03;     % "Resolution level" for computing surface coverage
        
        cylinder_start_matrix       = cylinder_centre_matrix - cylinder_length_list/2 .* cylinder_axis_matrix;
        [number_cylinders, num_dim] = size(cylinder_start_matrix);
        cylinder_surfcov_list       = zeros(number_cylinders, 1);

        for c = 1 : number_cylinders
            % Cylinder data
            cyl_point_cloud = cyl_point_cloud_cell{c};
            cylinder_axis   = cylinder_axis_matrix(c, :);
            cylinder_start  = cylinder_start_matrix(c, :);
            cylinder_radius = cylinder_radius_list(c);
            cylinder_length = cylinder_length_list(c);

            % Surface coverage
            number_points = size(cyl_point_cloud, 1);

            if number_points <= 1
                surfcov = 0;        % Surface coverage can't be calculated for fewer than two points
            else
                number_layers       = max(3, ceil(cylinder_length/res));
                number_sectors      = ceil(2*pi*cylinder_radius/res);
                number_sectors      = min(36, max(number_sectors, 8));
                min_distance        = 0.8 * cylinder_radius;
                [surfcov, ~, ~, ~]  = surface_coverage(cyl_point_cloud, cylinder_axis, cylinder_start, number_layers, number_sectors, min_distance);
            end

            cylinder_surfcov_list(c) = surfcov;
        end
 
    %% Radius correction %%
        % Unique branches
        branch_index_list   = unique(cylinder_branch_list);
        number_branches     = length(branch_index_list);

        % Update the geometry
        corr_cyl_radius_list        = zeros(number_cylinders, 1);
        corr_cyl_length_list        = zeros(number_cylinders, 1);
        corr_cyl_start_matrix       = zeros(number_cylinders, num_dim);
        corr_cyl_axis_matrix        = zeros(number_cylinders, num_dim);

        for n = 1 : number_branches
            %--% Branch cylinder data %--%
            % This branch's cylinders in order
            branch_index        = branch_index_list(n);
            branch_cylinders    = find(cylinder_branch_list == branch_index);

            branch_cyl_pos_list = cylinder_position_list(branch_cylinders);
            [~, order]          = sort(branch_cyl_pos_list);                       
            branch_cylinders    = branch_cylinders(order);

            % Their point clouds
            branch_cyl_point_cloud_cell = cyl_point_cloud_cell(branch_cylinders);

            % Geometry structure
            Branch_Cylinders = struct('radius', cylinder_radius_list(branch_cylinders), 'radius0', cylinder_radius_list(branch_cylinders), 'length', cylinder_length_list(branch_cylinders), 'start', cylinder_start_matrix(branch_cylinders, :), 'axis', cylinder_axis_matrix(branch_cylinders, :), 'BranchOrder', cylinder_BO_list(branch_cylinders), 'SurfCov', cylinder_surfcov_list(branch_cylinders));
            
            %--% Parent cylinder %--%
            first_cylinder_index    = branch_cylinders(1);
            branch_cyl_BO           = cylinder_BO_list(first_cylinder_index);
            branch_cyl_parent       = cylinder_parent_list(first_cylinder_index);

            if branch_cyl_BO == 0
                Parent_Cylinder = struct('radius', []);         % No parent cylinder for the stem
            else
                Parent_Cylinder = struct('radius', cylinder_radius_list(branch_cyl_parent), 'start', cylinder_start_matrix(branch_cyl_parent, :), 'axis', cylinder_axis_matrix(branch_cyl_parent, :), 'length', cylinder_length_list(branch_cyl_parent));
            end

            %--% Branch correction %--%
            Corr_Branch_Cylinders = adjustments(Branch_Cylinders, Parent_Cylinder, TreeQSM_Inputs, branch_cyl_point_cloud_cell);

            corr_cyl_radius_list(branch_cylinders)      = Corr_Branch_Cylinders.radius;
            corr_cyl_length_list(branch_cylinders)      = Corr_Branch_Cylinders.length;
            corr_cyl_start_matrix(branch_cylinders, :)  = Corr_Branch_Cylinders.start;
            corr_cyl_axis_matrix(branch_cylinders, :)   = Corr_Branch_Cylinders.axis;

            %--% Diagnostics plot %--%
            if Diagnostics == true
                % Cylinder colour map
                number_branch_cylinders = length(branch_cylinders);
                cylinder_cmap           = cbrewer('qual', 'Set1', max(number_branch_cylinders, 3));
                cylinder_cmap           = max(cylinder_cmap, 0);
                cylinder_cmap           = min(cylinder_cmap, 1);

                % Coordinates per cylinder
                number_coord = 1e2;

                % Cell array for the old and corrected branch cylinders
                Branch_Cylinders_cell   = {Branch_Cylinders, Corr_Branch_Cylinders};
                labels_cell             = {'Old', 'Corrected'};
                number_sets             = length(labels_cell);

                %--% Plot %--%
                figure(1)
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])     
    
                % Branch sets
                for s = 1 : number_sets
                    % This set's data
                    Set_Cylinders   = Branch_Cylinders_cell{s};
                    set_label       = labels_cell{s};

                    % Its subplot
                    subplot(1, number_sets, s);
                    hold on
                    grid on

                    title(set_label);

                    % Branch cylinders
                    for c = 1 : number_branch_cylinders
                        % This cylinder's data
                        cylinder_colour = cylinder_cmap(c, :);
                        cylinder_radius = Set_Cylinders.radius(c);
                        cylinder_length = Set_Cylinders.length(c);
                        cylinder_axis   = Set_Cylinders.axis(c, :);
                        cylinder_start  = Set_Cylinders.start(c, :);
                        cylinder_centre = cylinder_start + cylinder_length/2 * cylinder_axis;

                        % Its surface
                        [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, cylinder_length, cylinder_centre, cylinder_axis, number_coord);
                        surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', cylinder_colour, 'FaceAlpha', 0.25, 'LineWidth', 2);
                    end

                    % Parent cylinder
                    if branch_cyl_BO ~= 0
                        % Parent cylinder's geometry
                        cylinder_radius = cylinder_radius_list(branch_cyl_parent);
                        cylinder_length = cylinder_length_list(branch_cyl_parent);
                        cylinder_start  = cylinder_start_matrix(branch_cyl_parent, :);
                        cylinder_axis   = cylinder_axis_matrix(branch_cyl_parent, :);
                        cylinder_centre = cylinder_start + cylinder_length/2 * cylinder_axis;

                        % Its surface
                        [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, cylinder_length, cylinder_centre, cylinder_axis, number_coord);
                        surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.25, 'LineWidth', 2);
                    end

                    % Axes
                    xlabel('x [m]')
                    ylabel('y [m]')
                    zlabel('z [m]')
                
                    axis equal
                    view(45, 45);
                            
                    % Font size
                    set(gca, 'FontSize', 15);
                    set(gca, 'LineWidth', 2);
                
                    hold off
                end

                % Pause message
                disp('The branch cylinders were corrected. The script continues and figure closes upon a key-press.');
                pause();

                close(1);
            end
        end

        % NaN entries are not replaced
        NaN_bool                        = isnan(corr_cyl_radius_list);
        corr_cyl_radius_list(NaN_bool)  = cylinder_radius_list(NaN_bool);

        % Update the structure
        corr_cyl_centre_matrix  = corr_cyl_start_matrix + corr_cyl_length_list/2 .* corr_cyl_axis_matrix;
        corr_cyl_volume_list    = pi*corr_cyl_radius_list.^2 .* corr_cyl_length_list;

        Cylinders.UnmodRadius   = cylinder_radius_list;       % Like TreeQSM, the old radius is saved
        Cylinders.radius        = corr_cyl_radius_list;
        Cylinders.length        = corr_cyl_length_list;
        Cylinders.start         = corr_cyl_start_matrix;
        Cylinders.centre        = corr_cyl_centre_matrix;
        Cylinders.axis          = corr_cyl_axis_matrix;
        Cylinders.volume        = corr_cyl_volume_list;

    %% Local functions %%
        % Radius adjustment function from TreeQSM
        function Branch_Cylinders = adjustments(Branch_Cylinders, Parent_Cylinder, TreeQSM_Inputs, branch_cyl_point_cloud_cell)

            nc = size(Branch_Cylinders.radius,1);
            Mod = false(nc,1); % cylinders modified
            SC = Branch_Cylinders.SurfCov;
            
            %% Determine the maximum and the minimum radius
            % The maximum based on parent branch
            if ~isempty(Parent_Cylinder.radius)
              MaxR = 0.95*Parent_Cylinder.radius;
              MaxR = max(MaxR,TreeQSM_Inputs.MinCylRad);
            else
              % use the maximum from the bottom cylinders
              a = min(3,nc);
              MaxR = 1.25*max(Branch_Cylinders.radius(1:a));
            end
            
            MinR = min(Branch_Cylinders.radius(SC > 0.7));
            if ~isempty(MinR) && min(Branch_Cylinders.radius) < MinR/2
              MinR = min(MinR,min(Branch_Cylinders.radius(SC > 0.4)));
            elseif isempty(MinR)
              MinR = TreeQSM_Inputs.MinCylRad;
            end
            
            %% Check maximum and minimum radii
            I = Branch_Cylinders.radius < MinR;
            Branch_Cylinders.radius(I) = MinR;
            Mod(I) = true;
            if TreeQSM_Inputs.ParentCor || nc <= 3
              I = (Branch_Cylinders.radius > MaxR & SC < 0.7) | (Branch_Cylinders.radius > 1.2*MaxR);
              Branch_Cylinders.radius(I) = MaxR;
              Mod(I) = true;
              % For short branches modify with more restrictions
              if nc <= 3
                I = (Branch_Cylinders.radius > 0.75*MaxR & SC < 0.7);
                if any(I)
                  r = max(SC(I)/0.7.*Branch_Cylinders.radius(I),MinR);
                  Branch_Cylinders.radius(I) = r;
                  Mod(I) = true;
                end
              end
            end
            
            %% Use taper correction to modify radius of too small and large cylinders
            % Adjust radii if a small SurfCov and high SurfCov in the previous and
            % following cylinders
            for i = 2:nc-1
              if SC(i) < 0.7 && SC(i-1) >= 0.7 && SC(i+1) >= 0.7
                Branch_Cylinders.radius(i) = 0.5*(Branch_Cylinders.radius(i-1)+Branch_Cylinders.radius(i+1));
                Mod(i) = true;
              end
            end
            
            %% Use taper correction to modify radius of too small and large cylinders
            if TreeQSM_Inputs.TaperCor
              if max(Branch_Cylinders.radius) < 0.001
            
                %% Adjust radii of thin branches to be linearly decreasing
                if nc > 2
                  r = sort(Branch_Cylinders.radius);
                  r = r(2:end-1);
                  a = 2*mean(r);
                  if a > max(r)
                    a = min(0.01,max(r));
                  end
                  b = min(0.5*min(Branch_Cylinders.radius),0.001);
                  Branch_Cylinders.radius = linspace(a,b,nc)';
                elseif nc > 1
                  r = max(Branch_Cylinders.radius);
                  Branch_Cylinders.radius = [r; 0.5*r];
                end
                Mod = true(nc,1);
            
              elseif nc > 4
                %% Parabola adjustment of maximum and minimum
                % Define parabola taper shape as maximum (and minimum) radii for
                % the cylinders with low surface coverage
                branchlen = sum(Branch_Cylinders.length(1:nc)); % branch length
                L = Branch_Cylinders.length/2+[0; cumsum(Branch_Cylinders.length(1:nc-1))];
                Taper = [L; branchlen];
                Taper(:,2) = [1.05*Branch_Cylinders.radius; MinR];
                sc = [SC; 1];
            
                % Least square fitting of parabola to "Taper":
                A = [sum(sc.*Taper(:,1).^4) sum(sc.*Taper(:,1).^2); ...
                  sum(sc.*Taper(:,1).^2) sum(sc)];
                y = [sum(sc.*Taper(:,2).*Taper(:,1).^2); sum(sc.*Taper(:,2))];
                warning off
                x = A\y;
                warning on
                x(1) = min(x(1),-0.0001); % tapering from the base to the tip
                Ru = x(1)*L.^2+x(2); % upper bound parabola
                Ru( Ru < MinR ) = MinR;
                if max(Ru) > MaxR
                  a = max(Ru);
                  Ru = MaxR/a*Ru;
                end
                Rl = 0.75*Ru; % lower bound parabola
                Rl( Rl < MinR ) = MinR;
            
                % Modify radii based on parabola:
                % change values larger than the parabola-values when SC < 70%:
                I = Branch_Cylinders.radius > Ru & SC < 0.7;
                Branch_Cylinders.radius(I) = Ru(I)+(Branch_Cylinders.radius(I)-Ru(I)).*SC(I)/0.7;
                Mod(I) = true;
                % change values larger than the parabola-values when SC > 70% and
                % radius is over 33% larger than the parabola-value:
                I = Branch_Cylinders.radius > 1.333*Ru & SC >= 0.7;
                Branch_Cylinders.radius(I) = Ru(I)+(Branch_Cylinders.radius(I)-Ru(I)).*SC(I);
                Mod(I) = true;
                % change values smaller than the downscaled parabola-values:
                I = (Branch_Cylinders.radius < Rl & SC < 0.7) | (Branch_Cylinders.radius < 0.5*Rl);
                Branch_Cylinders.radius(I) = Rl(I);
                Mod(I) = true;
            
              else
                %% Adjust radii of short branches to be linearly decreasing
                R = Branch_Cylinders.radius;
                if nnz(SC >= 0.7) > 1
                  a = max(R(SC >= 0.7));
                  b = min(R(SC >= 0.7));
                elseif nnz(SC >= 0.7) == 1
                  a = max(R(SC >= 0.7));
                  b = min(R);
                else
                  a = sum(R.*SC/sum(SC));
                  b = min(R);
                end
                Ru = linspace(a,b,nc)';
                I = SC < 0.7 & ~Mod;
                Branch_Cylinders.radius(I) = Ru(I)+(R(I)-Ru(I)).*SC(I)/0.7;
                Mod(I) = true;
            
              end
            end
            
            %% Modify starting points by optimising them for given radius and axis
            nr = size(branch_cyl_point_cloud_cell,1);
            for i = 1:nc
              if Mod(i)
                if nr == nc
                  Reg = branch_cyl_point_cloud_cell{i};

                  if isempty(Reg)
                      continue
                  end
                elseif i > 1
                  Reg = branch_cyl_point_cloud_cell{i-1};

                  if isempty(Reg)
                      continue
                  end
                end

                if abs(Branch_Cylinders.radius(i)-Branch_Cylinders.radius0(i)) > 0.005 && ...
                    (nr == nc || (nr < nc && i > 1))
                  P = Reg-Branch_Cylinders.start(i,:);
                  [U,V] = orthonormal_vectors(Branch_Cylinders.axis(i,:));
                  P = P*[U V];
                  cir = least_squares_circle_centre(P,[0 0],Branch_Cylinders.radius(i));
                  if cir.conv && cir.rel
                    Branch_Cylinders.start(i,:) = Branch_Cylinders.start(i,:)+cir.point(1)*U'+cir.point(2)*V';
                    Branch_Cylinders.mad(i,1) = cir.mad;
                    [~,V,h] = distances_to_line(Reg,Branch_Cylinders.axis(i,:),Branch_Cylinders.start(i,:));
                    if min(h) < -0.001
                      Branch_Cylinders.length(i) = max(h)-min(h);
                      Branch_Cylinders.start(i,:) = Branch_Cylinders.start(i,:)+min(h)*Branch_Cylinders.axis(i,:);
                      [~,V,h] = distances_to_line(Reg,Branch_Cylinders.axis(i,:),Branch_Cylinders.start(i,:));
                    end
                    a = max(0.02,0.2*Branch_Cylinders.radius(i));
                    nl = ceil(Branch_Cylinders.length(i)/a);
                    nl = max(nl,4);
                    ns = ceil(2*pi*Branch_Cylinders.radius(i)/a);
                    ns = max(ns,10);
                    ns = min(ns,36);
                    Branch_Cylinders.SurfCov(i,1) = surface_coverage2(...
                      Branch_Cylinders.axis(i,:),Branch_Cylinders.length(i),V,h,nl,ns);
                  end
                end
              end
            end
            
            %% Continuous branches
            % Make cylinders properly "continuous" by moving the starting points
            % Move the starting point to the plane defined by parent cylinder's top
            % if nc > 1
            %   for j = 2:nc
            %     U = Branch_Cylinders.start(j,:)-Branch_Cylinders.start(j-1,:)-Branch_Cylinders.length(j-1)*Branch_Cylinders.axis(j-1,:);
            %     if (norm(U) > 0.0001)
            %       % First define vector V and W which are orthogonal to the
            %       % cylinder axis N
            %       N = Branch_Cylinders.axis(j,:)';
            %       if norm(N) > 0
            %         [V,W] = orthonormal_vectors(N);
            %         % Now define the new starting point
            %         x = [N V W]\U';
            %         Branch_Cylinders.start(j,:) = Branch_Cylinders.start(j,:)-x(1)*N';
            %         if x(1) > 0
            %           Branch_Cylinders.length(j) = Branch_Cylinders.length(j)+x(1);
            %         elseif Branch_Cylinders.length(j)+x(1) > 0
            %           Branch_Cylinders.length(j) = Branch_Cylinders.length(j)+x(1);
            %         end
            %       end
            %     end
            %   end
            % end
            
            %% Connect far away first cylinder to the parent
            if ~isempty(Parent_Cylinder.radius)
              [d,V,h,B] = distances_to_line(Branch_Cylinders.start(1,:),Parent_Cylinder.axis,Parent_Cylinder.start);
              d = d-Parent_Cylinder.radius;
              if d > 0.001
                taper = Branch_Cylinders.start(1,:);
                E = taper+Branch_Cylinders.length(1)*Branch_Cylinders.axis(1,:);
                V = Parent_Cylinder.radius*V/norm(V);
                if h >= 0 && h <= Parent_Cylinder.length
                  Branch_Cylinders.start(1,:) = Parent_Cylinder.start+B+V;
                elseif h < 0
                  Branch_Cylinders.start(1,:) = Parent_Cylinder.start+V;
                else
                  Branch_Cylinders.start(1,:) = Parent_Cylinder.start+Parent_Cylinder.length*Parent_Cylinder.axis+V;
                end
                Branch_Cylinders.axis(1,:) = E-Branch_Cylinders.start(1,:);
                Branch_Cylinders.length(1) = norm(Branch_Cylinders.axis(1,:));
                Branch_Cylinders.axis(1,:) = Branch_Cylinders.axis(1,:)/Branch_Cylinders.length(1);
              end
            end
            end

end

