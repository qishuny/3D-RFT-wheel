%%
classdef HmapCellContext < handle
    
    properties
        HmapCells;
        activeQueue = Queue;
        minx; miny;
        maxx; maxy;
        resolution; 
        maxheightchange_laststep;
        dt = 0.1;
    end
    properties (Constant)
        flow_rate = 0.01;
        angle_repose = 0.6;
        mycolormap = [
            0.5159    0.3224    0.2053  %my colormap for visualization
            0.5357    0.3348    0.2132
            0.5556    0.3472    0.2211
            0.5754    0.3596    0.2290
            0.5952    0.3720    0.2369
            0.6151    0.3844    0.2448
            0.6349    0.3968    0.2527
            0.6548    0.4092    0.2606
            0.6746    0.4216    0.2685
            0.6944    0.4340    0.2764
            0.7143    0.4464    0.2843
            0.7341    0.4588    0.2922
            0.7540    0.4712    0.3001
            0.7738    0.4836    0.3080
            0.7937    0.4960    0.3159
            0.8135    0.5084    0.3238
            0.8333    0.5208    0.3317
            0.8532    0.5332    0.3396
            0.8730    0.5456    0.3475
            0.8929    0.5580    0.3554
            0.9127    0.5704    0.3633
            0.9325    0.5828    0.3712
            0.9524    0.5952    0.3790
            0.9722    0.6076    0.3869
            0.9921    0.6200    0.3948
            1.0000    0.6324    0.4027
            1.0000    0.6448    0.4106
            1.0000    0.6572    0.4185
            1.0000    0.6696    0.4264
            1.0000    0.6820    0.4343
            1.0000    0.6944    0.4422
            1.0000    0.7068    0.4501
            1.0000    0.7192    0.4580
            1.0000    0.7316    0.4659
            1.0000    0.7440    0.4738
            1.0000    0.7564    0.4817
            1.0000    0.7688    0.4896
            1.0000    0.7812    0.4975];
    end
    %%
    methods
        % constructor
        function this = HmapCellContext(minx,miny,maxx,maxy,resolution)
            if nargin == 0
                
            elseif nargin == 5
                this.resolution = resolution;
                this.minx = minx; this.miny = miny; 
                this.maxx = maxx; this.maxy = maxy;
                if (maxx-minx) ~= (maxy-miny)
                    disp('Careful! Your x size and y size are different')
                end
%                 cellxsize = (maxx-minx)/(resolution-1);
%                 cellysize = (maxy-miny)/(resolution-1);
                % empty cell array to hold my HmapCells 
                this.HmapCells = cell(resolution,resolution);
                % assign i,j,x,y and neighbors to all my HmapCells
                this.initialize();
            end  
        end
        
        function dx = get_dx(this)
            dx = (this.maxx - this.minx)/(this.resolution-1);
        end
        
        
        function getActivelist(this)
            for i = 1:size(this.HmapCells,1)
                for j = 1:size(this.HmapCells,2)
                    if this.HmapCells{i,j}.isActive
                        this.activeHmapCells{end+1} = this.HmapCells{i,j};
                    end
                end
            end            
        end
        
        
        function initialize(this)
            x_inc = (this.maxx-this.minx)/(this.resolution-1);
            y_inc = (this.maxy-this.miny)/(this.resolution-1);
            for i = 1:size(this.HmapCells,1)
                for j = 1:size(this.HmapCells,2)
                    this.HmapCells{i,j} = HmapCell(i,j,(i-1)*x_inc+this.minx,(j-1)*y_inc+this.miny);
                    this.HmapCells{i,j}.height = 0;
                end
            end
            
            %set back pointers for neighbors (4 connectedness)
            for i = 1:size(this.HmapCells,1)
                for j = 1:size(this.HmapCells,2)
                    if i>1
                        this.HmapCells{i,j}.neighbors{end+1} = this.HmapCells{i-1,j};
                    end
                    if j>1
                        this.HmapCells{i,j}.neighbors{end+1} = this.HmapCells{i  ,j-1};
                    end
                    if j<size(this.HmapCells,2)
                        this.HmapCells{i,j}.neighbors{end+1} = this.HmapCells{i  ,j+1};
                    end
                    if i<size(this.HmapCells,1)
                        this.HmapCells{i,j}.neighbors{end+1} = this.HmapCells{i+1,j};
                    end
                end
            end
        end
            
        
        
        function updateOnestepActivegrid(this,dt)
            if nargin == 1
                dt = 0.1;
            end
            dx = this.get_dx();
            k = this.flow_rate;
            angleRepose = this.angle_repose;
            
%             this.activeQueue.list
%             this.activeQueue
            % do update for each active grid
            for i = 1:length(this.activeQueue.list)
                
                h_cell = this.activeQueue.list{i}; %copy of reference to it
                sand_from_neighbors = 0; %this will contain flow from all neighbors for this h_grid
                
                avalancheTriggered = false;
                % sand flow q from all neighbors
                for neigh_idx = 1:length(h_cell.neighbors)
                    neighbor = h_cell.neighbors{neigh_idx};
                    slope = (neighbor.height - h_cell.height)/dx;
                    if abs(slope) > angleRepose
                        %if flow from any node this node is now active
                        avalancheTriggered = true;
                        % and add this neighbor to the active queue
                        this.activeQueue.add(h_cell.neighbors{neigh_idx});
                    end
                    sand_from_neighbors = sand_from_neighbors + ...
                        k*dx*sign(slope)*(abs(slope)-angleRepose)*(abs(slope)>angleRepose);
                end
                % then this node is active for the next round
                h_cell.isActive = avalancheTriggered;
                % sand flow rate q * dt = sand flow volume
                delta_h = sand_from_neighbors*dt/(dx*dx);
                % save new height on height_temp, update all at once later
                h_cell.height_temp = h_cell.height + delta_h;
            end
            
            [hmap, hmap_temp] = this.get_Matrix();
            this.maxheightchange_laststep = max(max(hmap_temp - hmap)); %save how much changed for flag
            % update all at once
            for i = 1:length(this.activeQueue.list)
                this.activeQueue.list{i}.height = this.activeQueue.list{i}.height_temp;
            end
        end
        
        

        function updateOnestepAllgrid(this,dt)
            % update all grid and add Active queue for the nodes where
            % avalanche happened and its neighbors run this at the very
            % beginning of first impulse
            this.activeQueue.clear(); %we will start with clean list of activeQueue
            if nargin == 1
                dt = 0.1;
            end
            dx = this.get_dx();
            k = this.flow_rate;
            angleRepose = this.angle_repose;
            
            for i = 1:size(this.HmapCells,1)
                for j = 1:size(this.HmapCells,2)
                    
                    h_cell = this.HmapCells{i,j}; %copy of it
                    sand_from_neighbors = 0; %this will contain flow from all neighbors for this h_grid
                    
                    avalancheTriggered = false;
                    % sand flow q from all neighbors
                    for neigh_idx = 1:length(h_cell.neighbors)
                        neighbor = h_cell.neighbors{neigh_idx};
                        slope = (neighbor.height - h_cell.height)/dx;
                        if abs(slope) > angleRepose
                            %if flow from any node this node is now active
                            avalancheTriggered = true;
                            % add the neighbor to the activeQueue
                            this.activeQueue.add(h_cell.neighbors{neigh_idx});
%                             disp('activeQueue added')
%                             h_cell.neighbors{neigh_idx}
%                             this.activeQueue.list{end}
                        end
                        sand_from_neighbors = sand_from_neighbors + ...
                            k*dx*sign(slope)*(abs(slope)-angleRepose)*(abs(slope)>angleRepose);
                    end
                    % then this node is active for the next round
                    h_cell.isActive = avalancheTriggered;
                    % sand flow rate q * dt = sand flow volume
                    delta_h = sand_from_neighbors*dt/(dx*dx);
                    % save new height on height_temp, update all at once later
                    h_cell.height_temp = h_cell.height + delta_h;
                end
            end
            
            [hmap, hmap_temp] = this.get_Matrix();
            this.maxheightchange_laststep = max(max(hmap_temp - hmap)); %save how much changed for flag
            % update all at once
            for i = 1:size(this.HmapCells,1)
                for j = 1:size(this.HmapCells,2)
                    this.HmapCells{i,j}.height = this.HmapCells{i,j}.height_temp;
                end
            end
        end
        
        function [hmatrix,htempmatrix] = get_Matrix(this)
            hmatrix = zeros(size(this.HmapCells));
            htempmatrix = zeros(size(hmatrix));
            for i = 1:size(this.HmapCells,1)
                for j = 1:size(this.HmapCells,2)
                    hmatrix(i,j) = this.HmapCells{i,j}.height;
                    try
                        htempmatrix(i,j) = this.HmapCells{i,j}.height_temp;
                    catch
                        htempmatrix(i,j) = 0; 
                    end
                end
            end
        end
        
        % update hmap until maximum volume flown last step is less then
        % threshold (volume_stop_thr)
        function steadyState(this,dt,doDynamicsVisualization)
            if nargin == 1 %i.e) no input
                dt = this.dt;
            end
            this.updateOnestepAllgrid(dt);
            maxchange = this.maxheightchange_laststep;
            volume_stop_thr = 0.000001;  %0.0000001
%             maxchange*this.get_dx()^2 > volume_stop_thr
            
            while maxchange*this.get_dx()^2 > volume_stop_thr
                this.updateOnestepActivegrid(dt);
%                 this.updateOnestepAllgrid(dt);
                maxchange = this.maxheightchange_laststep;
                maxchange*this.get_dx()^2
%                 activelist = length(this.activeQueue.list)
%                 vol_change = maxchange*this.get_dx()^2;
                if doDynamicsVisualization
                    VisualizeContext.smethodsVisualize(this);
                end
%                 maxchange = 100;
            end
            VisualizeContext.smethodsVisualize(this);
        end
        
        
        function resolve_collision(this,b_idx,b_z,v)
            
            % b_idx: blade index
            % b_z : corresponding height 
            % v : velocity of the blade v = [vx; vy];
            if v ~= 0
                v = v/(sqrt(v(1)^2+v(2)^2)); %make it into unit vector
            end
            this.blocked_idx = b_idx;
            while collision_detected(this.matrix,b_idx,b_z) %check if there is collision
                % if there is resolve it
                for i = 1:size(b_idx,2)
                    collision_depth = b_z(i) - this.matrix(b_idx(1,i),b_idx(2,i));
                    if collision_depth < 0 % if there is collision
                        % remove collision amount of sand
                        this.matrix(b_idx(1,i),b_idx(2,i)) = b_z(i);
                        % and move that sand to the direction of v
                        if v == 0 % if it is not moving, move it to random direction
                            th = rand;
                            v(1) = cos(th);  v(2) = sin(th);
                            d_idx(1) = round(b_idx(1,i) + v(1));
                            d_idx(2) = round(b_idx(2,i) + v(2)); %displaced index
                        else
                            d_idx(1) = round(b_idx(1,i) + v(1));
                            d_idx(2) = round(b_idx(2,i) + v(2)); %displaced index
                        end
                        this.matrix(d_idx(1),d_idx(2)) = this.matrix(d_idx(1),d_idx(2)) - collision_depth;
                        
                    end
                end
            end
        end
        
        
    end
end