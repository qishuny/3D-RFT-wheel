classdef Heightmap < handle
    %%
    properties
        n; %grids per meter
        dx;
        matrix;
        min_corner;
        max_corner;
        blade_idx;
        blade_idx_front;
        blade_idx_back;
        dT_default; % default step time for visualization
        mycolormap;
    end
    
    methods
        

        % if there is no specification we go with this
        function this = Heightmap()
            this.n = 10;
            this.dx = 1/this.n;
            this.min_corner = [0 0]';
            this.max_corner = [4 4]';
            this.matrix = zeros(10*4+1);
            this.dT_default = 0.1;
            this.mycolormap = [
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
%             this.blade_idx = this.get_coord(this,[1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8;...
%                                                   1 1   1   1   1   1   1   1   1  ]);
        end

        function set_param(this,n,minx,maxx)
            this.n = n;
            this.dx = 1/n;
            this.min_corner = [minx minx];
            this.max_corner = [maxx maxx];
            this.matrix = zeros(int64((maxx-minx)/this.dx)+1);           
        end
        

        function synchronize(this,blade_pos,blade_pos_front,blade_pos_back,blade_z)
        % synchronize with the blade position
        % block flow where the blade and blade back is located
        % remove sand from blade and add sand to blade pos front
        % amount of sand removed and added depends on the depth
            if nargin <4
                blade_z = zeros(size(blade_pos));
            end
            if length(blade_z) < 2
                % if blade_z is just a scalar, make it into vector for
                % later use
                blade_z = blade_z * ones(size(blade_pos));
            end
            this.blade_idx = this.get_coord(this,blade_pos);
            this.blade_idx_front = this.get_coord(this,blade_pos_front);
            this.blade_idx_back = this.get_coord(this,blade_pos_back);
            for i = 1:length(this.blade_idx)
                x_idx = this.blade_idx(1,i);
                y_idx = this.blade_idx(2,i);
                x_idx_f = this.blade_idx_front(1,i);
                y_idx_f = this.blade_idx_front(2,i);
                this.matrix(x_idx_f,y_idx_f) = this.matrix(x_idx_f,y_idx_f) + (this.matrix(x_idx,y_idx) - blade_z(i)); % z is negative
                % add amount of sand displaced to blade_f
                this.matrix(x_idx,y_idx) = blade_z(i);
                % remove amound of sand displaced in blade
            end
            
        end
        
%         function push_sand(this)
%             for i = 1:length(this.blade_idx)
%                 x_idx = this.blade_idx(1,i);
%                 y_idx = this.blade_idx(2,i);
%                 x_idx_f = this.blade_idx_front(1,i);
%                 y_idx_f = this.blade_idx_front(2,i);
%                 
%                 this.matrix(x_idx_f,y_idx_f) = this.matrix(x_idx_f,y_idx_f) + this.matrix(x_idx,y_idx);
%                 this.matrix(x_idx,y_idx) = 0;
%             end            
%         end
    
        function add_sand(this,x,y,height)
            coord = this.get_coord(this,[x;y]);
            this.matrix(coord(1),coord(2)) = this.matrix(coord(1),coord(2)) + height;
        end
        
        function steady_state(this)

%             maxchange = max(max(this.matrix - mtx_method_update_height(this.matrix,this.dx,temp,this.dT_default)));          
            temp = [this.blade_idx this.blade_idx_back]; %1st row x 2nd row y
            maxchange = max(max(this.matrix - mtx_method_update_height(this.matrix,this.dx,temp,this.dT_default)));          
            volume_stop_thr = 0.0000001;
            while maxchange*this.dx^2 > volume_stop_thr
                matrix_old = this.matrix;
                temp = [this.blade_idx this.blade_idx_back];
                this.matrix = mtx_method_update_height(this.matrix,this.dx,temp,this.dT_default);
                maxchange = max(max(matrix_old - this.matrix));
            end

%             for i = 1:200
%                 temp = [this.blade_idx this.blade_idx_back];
%                 this.matrix = mtx_method_update_height(this.matrix,this.dx,temp,this.dT_default);
%             end
%             while maxchange*this.dx^2 > 0.0000001
%                 matrix_old = this.matrix;
%                 temp = [this.blade_idx this.blade_idx_back];
%                 this.matrix = mtx_method_update_height(this.matrix,this.dx,temp,this.dT_default);
%                 maxchange = max(max(matrix_old - this.matrix));
%             end
        end
        
        function update(this,steps)
           temp = [this.blade_idx this.blade_idx_back];
           for i = 1:steps
                temp = [this.blade_idx this.blade_idx_back];
                this.matrix = mtx_method_update_height(this.matrix,this.dx,temp,this.dT_default);
           end         
        end
        
        function visualize(this)
            terrain_sizex = this.max_corner(1)-this.min_corner(1);
            xmin = this.min_corner(1);
            xmax = this.max_corner(1);
            gridx = linspace(xmin,xmax,this.n*terrain_sizex+1)';
            gridX = [];
            for i = 1:terrain_sizex*this.n+1
                gridX = [gridX gridx];
            end
            
            terrain_sizey = this.max_corner(2)-this.min_corner(2);
            ymin = this.min_corner(2);
            ymax = this.max_corner(2);
            gridy = linspace(ymin,ymax,this.n*terrain_sizey+1);
            gridY = [];
            for i = 1:terrain_sizey*this.n+1
                gridY = [gridY; gridy];
            end

            temp = [this.blade_idx this.blade_idx_back];
            this.matrix = mtx_method_update_height(this.matrix,this.dx,temp,this.dT_default);
            s = surf(gridX,gridY,this.matrix);
            s.EdgeColor = 'none';
            colormap(this.mycolormap);
            axis([xmin xmax ymin ymax -xmax/10 xmax])
            view([290 30])
            drawnow
        end
        
        function visualize_dynamics(this)
            
            terrain_sizex = this.max_corner(1)-this.min_corner(1);
            xmin = this.min_corner(1);
            xmax = this.max_corner(1);
            gridx = linspace(xmin,xmax,this.n*terrain_sizex+1)';
            gridX = [];
            for i = 1:terrain_sizex*this.n+1
                gridX = [gridX gridx];
            end
            
            terrain_sizey = this.max_corner(2)-this.min_corner(2);
            ymin = this.min_corner(2);
            ymax = this.max_corner(2);
            gridy = linspace(ymin,ymax,this.n*terrain_sizey+1);
            gridY = [];
            for i = 1:terrain_sizey*this.n+1
                gridY = [gridY; gridy];
            end
            
            temp = [this.blade_idx this.blade_idx_back];
            maxchange = max(max(this.matrix - mtx_method_update_height(this.matrix,this.dx,temp,this.dT_default)));          
            rep_i =  0;
            volume_stop_thr = 0.0000001;
            while maxchange*this.dx^2 > volume_stop_thr
                rep_i = rep_i + 1;
                matrix_old = this.matrix;
                temp = [this.blade_idx this.blade_idx_back];
                this.matrix = mtx_method_update_height(this.matrix,this.dx,temp,this.dT_default);
                s = surf(gridX,gridY,this.matrix);
                s.EdgeColor = 'none';
                colormap(this.mycolormap);
                axis([xmin xmax ymin ymax -xmax/10 xmax])
                view([290 30])
                drawnow
                maxchange = max(max(matrix_old - this.matrix));
            end
%             rep_i

%             for i = 1:200
%                 temp = [this.blade_idx this.blade_idx_back];
%                 matrix_old = this.matrix;
%                 this.matrix = mtx_method_update_height(this.matrix,this.dx,temp,this.dT_default);
% %                 size(gridX)
% %                 size(gridY)
% %                 size(this.matrix)
%                 s = surf(gridX,gridY,this.matrix);
%                 s.EdgeColor = 'none';
%                 colormap(this.mycolormap);
%                 axis([xmin xmax ymin ymax -xmax/10 xmax])
%                 drawnow
% %                 refreshdata
%             end

%             [temp, x_indices] = max(matrix_old - this.matrix);
%             [max_update, y_idx] = max(temp);
%             x_idx = x_indices(y_idx);
%             [x_idx, y_idx vpa(max_update)]

        end
 
        function visualize_onestep(this,dT)
            
            [gridX, gridY] = make_grid(this.min_corner(1),this.max_corner(1),this.min_corner(2),this.max_corner(2),this.n);
            
            for i = 1:1
%                 change = this.matrix - mtx_method_update_height(this.matrix,this.dx);
%                 max(max(abs(change)))
%                 this.matrix = mtx_method_update_height(this.matrix,this.dx);
                temp = [this.blade_idx this.blade_idx_back];
                this.matrix = mtx_method_update_height(this.matrix,this.dx,temp,dT);
                s = surf(gridX,gridY,this.matrix);
                s.EdgeColor = 'none';
                colormap(this.mycolormap);
                axis([this.min_corner(1) this.max_corner(1) ...
                      this.min_corner(2) this.max_corner(2) ...
                      -this.min_corner(1)/10 this.max_corner(1)])
                view([290 30])
                drawnow
            end
        end        
        
        
        
    end
    
    methods(Static)
        function indices = get_coord(this,p)
            x_idx = int64(p(1,:)/this.dx)+1;
            y_idx = int64(p(2,:)/this.dx)+1;
            indices = [x_idx; y_idx];
        end
        
    end
    
    
end