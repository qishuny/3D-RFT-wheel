classdef Hmap < handle
    %%
    properties (Constant)
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
    properties
        n; %grids per meter
        matrix;
        minx; miny; maxx; maxy;
        dT;
        volume_stop_thr = 0.00001;
        reachedSteadystate = false;
        % from Vehicle
        blocked_idx;
        displaced_idx;
        optOn = false; 
        optIdx = [0 0; 0 0]; % xmin xmax; ymin ymax of the matrix indices I update
    end
    
    methods
        function this = Hmap(minx,miny,maxx,maxy,n)
            if nargin == 0
                this.n = 10;
                this.minx = 0; this.miny = 0; 
                this.maxx = 4; this.maxy = 4;
                M = round(this.n*(this.maxx-this.minx)+1);
                N = round(this.n*(this.maxy-this.miny)+1);
                this.matrix = zeros(M,N);
                this.dT = 2;
                this.optIdx = [1 M;1 N];
            elseif nargin == 5
                this.n = n;
                this.minx = minx; this.miny = miny; 
                this.maxx = maxx; this.maxy = maxy;
                M = round(this.n*(this.maxx-this.minx)+1);
                N = round(this.n*(this.maxy-this.miny)+1);
                this.matrix = zeros(M,N);
                this.dT = 2;
            else 
                error('input must be appropriate');
            end
        end
            
        function dx = get_dx(this)
            dx = 1/this.n;
        end      
        
        function update_onestep(this) 
            if this.optOn
                collision_idx = this.blocked_idx; 
                mapSize = size(this.matrix); 
                % additional area to cover during simulation
                optBorder = round(mapSize / 8); 
                % minx
                this.optIdx(1,1) = min(collision_idx(1,:)) - optBorder(1); 
                if this.optIdx(1,1) < 1
                    this.optIdx(1,1) = 1;
                end
                % miny
                this.optIdx(2,1) = min(collision_idx(2,:)) - optBorder(2); 
                if this.optIdx(2,1) < 1
                    this.optIdx(2,1) = 1;
                end
                % maxx
                this.optIdx(1,2) = max(collision_idx(1,:)) + optBorder(1); 
                if this.optIdx(1,2) > mapSize(1)
                    this.optIdx(1,2) = mapSize(1);
                end
                % maxy
                this.optIdx(2,2) = max(collision_idx(2,:)) + optBorder(2); 
                if this.optIdx(2,2) > mapSize(2)
                    this.optIdx(2,2) = mapSize(2);  
                end                
                xmin = this.optIdx(1,1); xmax = this.optIdx(1,2);
                ymin = this.optIdx(2,1); ymax = this.optIdx(2,2); 
                % new index for the corresponding matrix!
                bl_idx = this.blocked_idx - [xmin-1; ymin-1];
                
                this.matrix(xmin:xmax,ymin:ymax) = ...
                    mtx_method_update_height(this.matrix(xmin:xmax,ymin:ymax),this.get_dx,bl_idx,this.dT);
            else 
                matrix_old = this.matrix;
                this.matrix = mtx_method_update_height(this.matrix,this.get_dx,this.blocked_idx,this.dT);
                maxchange = max(max(matrix_old - this.matrix));
                if maxchange*this.get_dx^2 < this.volume_stop_thr
                    this.reachedSteadystate = true;
                else
                    this.reachedSteadystate = false;
                end     
            end
        end
        
        % need to upgrade, can't check collision everytime
        function Steadystate(this,collision_idx,collision_z,v)
            % preprocess for optimization 
            if this.optOn
                mapSize = size(this.matrix); 
                % additional area to cover during simulation
                optBorder = round(mapSize / 8); 
                % minx
                this.optIdx(1,1) = min(collision_idx(1,:)) - optBorder(1); 
                if this.optIdx(1,1) < 1
                    this.optIdx(1,1) = 1;
                end
                % miny
                this.optIdx(2,1) = min(collision_idx(2,:)) - optBorder(2); 
                if this.optIdx(2,1) < 1
                    this.optIdx(2,1) = 1;
                end
                % maxx
                this.optIdx(1,2) = max(collision_idx(1,:)) + optBorder(1); 
                if this.optIdx(1,2) > mapSize(1)
                    this.optIdx(1,2) = mapSize(1);
                end
                % maxy
                this.optIdx(2,2) = max(collision_idx(2,:)) + optBorder(2); 
                if this.optIdx(2,2) > mapSize(2)
                    this.optIdx(2,2) = mapSize(2);  
                end
            end
            
            this.resolve_collision(collision_idx,collision_z,v);
%             this.update_onestep;
%             this.resolve_collision(collision_idx,collision_z,v);
            
%             for steadystate=1:round(1/this.get_dx^2)
%                 this.update_onestep;
%                 this.resolve_collision(collision_idx,collision_z,v);
%             end
            for steadystate=1:200 % how do I pick a good number?
                this.update_onestep;
%                 this.resolve_collision(collision_idx,collision_z,v);
            end
        end        
 
%         function add_sand(this, xidx, yidx, height)
%             this.matrix(xidx, ydix) = this.matrix(xidx, yidx) + height; 
%         end
        
        % deprecate, get_coord is handled by IndexHandler
        function add_sand(this,x,y,height)
            ind_array = get_coord([x;y],this.get_dx,this.minx,this.miny);
            this.matrix(ind_array(1,:),ind_array(2,:)) = this.matrix(ind_array(1,:),ind_array(2,:)) + height;
        end
        
        function resolve_collision(this,collision_idx,collision_z,v)
            % collision_idx: blocked index
            % collision_z : corresponding height 
            % v : velocity of the blade v = [vx; vy];
            if norm(v) ~= 0
                v = v/(sqrt(v(1)^2+v(2)^2)); %make it into unit vector
            end
            % our new blocked idx is this collision idx
            this.blocked_idx = collision_idx;

            while collision_detected(this.matrix,collision_idx,collision_z) %check if there is collision
                % if there is resolve it
                for i = 1:size(collision_idx,2)
                    collision_depth = collision_z(i) - this.matrix(collision_idx(1,i),collision_idx(2,i));
                    if collision_depth < 0 % if there is collision
                        % if you want to make sure location is blocked only
                        % after collision happens (ex. back of an wheel)
                        % you have to add to the index
%                         this.blocked_idx = [this.blocked_idx collision_idx(:,i)];
                        % remove collision amount of sand
                        this.matrix(collision_idx(1,i),collision_idx(2,i)) = collision_z(i);
                        % and move that sand to the direction of v
                        if v == 0 % if it is not moving, move it to random direction
                            th = rand*2*pi;
                            vtemp(1) = cos(th);  vtemp(2) = sin(th);
                            d_idx(1) = round(collision_idx(1,i) + vtemp(1));
                            d_idx(2) = round(collision_idx(2,i) + vtemp(2)); %displaced index
                        else
%                             sprintf('v1 is %d, v2 is %d',v(1),v(2));
                            d_idx(1) = round(collision_idx(1,i) + v(1));
                            d_idx(2) = round(collision_idx(2,i) + v(2)); %displaced index
                        end
                        this.matrix(d_idx(1),d_idx(2)) = this.matrix(d_idx(1),d_idx(2)) - collision_depth;
                    end
                end
            end
        end       
    end 
    
    methods (Static)
        
    end
    
    
    
end

function detected = collision_detected(matrix,b_idx,b_z)
for i = 1:length(b_idx)
    collision_depth = b_z(i) - matrix(b_idx(1,i),b_idx(2,i));
    if collision_depth < 0
        detected = true;
        return
    end
end
detected = false;
end