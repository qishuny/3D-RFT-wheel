%% 

classdef Vehicle < handle
    properties
%         x; y; th; 
        Blade;
        BladeFront;
        HT;
        HT4;
        z_function;
        PosArray3; % need to implement this part
        PosArray3_displaced;
        PosArray3_cell; % cell of 
        
        PartsV={}; %cell of PartVisualization objects for visualization

        % from hmap
        dx; minx; miny; %for getting indices
    end

    methods 
        function this = Vehicle(x,y,th,dx,Blade,BladeFront)
            if nargin == 0
                this.dx = 0.1;
                this.Blade = Part(-0.3,-0.3,0.3,0.3);
                this.BladeFront = Part(-0.3,0.3,0.6,0.6);
                this.HT  = get_HT(0,0,0);
                this.HT4 = get_HT4(0,0,0,0);
                this.z = HeightFunction();
            elseif nargin == 4
                this.dx = dx;
                this.HT = get_HT(x,y,th);
                this.HT4 = get_HT4(x,y,0,th);
%                 this.Blade = Part(-0.3,-0.3,0.3,0.3);
%                 this.BladeFront = Part(-0.3,0.3,0.6,0.6);
            elseif nargin == 6
                this.dx = dx;
                this.HT = get_HT(x,y,th);
                this.HT4 = get_HT4(x,y,0,th);
                this.Blade = Blade;
                this.BladeFront = BladeFront;
            else 
                error('input must be appropriate');
            end
            this.PosArray3_cell = {};
            this.PartsV = {};
        end

        
        function x = get_x(this)
            x = this.HT(1,3);
        end
        function y = get_y(this)
            y = this.HT(2,3);
        end
        function th = get_th(this)
            th = atan2(this.HT(2,1),this.HT(1,1));
        end
        function invHT = get_invHT(this)
            th = this.get_th;
            x = this.get_x; y = this.get_y;
            R = [cos(th) -sin(th);
                sin(th) cos(th)];
            p = -R'*[x;y];
            invHT = [cos(th) sin(th) p(1)
                -sin(th) cos(th) p(2)
                0 0 1];
        end

        function synchronize(this,Hmap)
            this.dx = Hmap.get_dx;
            this.minx = Hmap.minx; 
            this.miny = Hmap.miny; 
        end
        
        function move(this,forward,turn)
            % forward is 1  backward is -1
            % turn 1 for ccw 10 deg, -1 for cw 10 deg
            t = turn*10*pi/180;
            p = forward*this.dx;
            R = [cos(t) -sin(t) 0;
                sin(t) cos(t)  p;
                0       0        1];
            R4 = [cos(t) -sin(t) 0 0;
                sin(t) cos(t)  0 p;
                0       0      1 0;
                0      0       0 1];
            this.HT = this.HT*R; 
            this.HT4 = this.HT4*R4;
        end
        
        function pos_array_cell = get_plot_pos(this)
            for i = 1:length(this.PartsV)
                for pos_array_idx = 1:length(this.PartsV{i}.PosArray)
                    pos_array_cell{i} =this.HT4 * this.PartsV{i}.pos_array;
                end
            end
        end
        
        function visualize(this)            
            hold on
            for i = 1:length(this.PartsV) %for wheel1, wheel2, wheel3, ...
                for idx_lines = 1:length(this.PartsV{i}.lines)%for each lines in wheel1
                    plot_points = this.HT4 * this.PartsV{i}.lines{idx_lines}.pos_array;
                    plot3(plot_points(1,:),plot_points(2,:),plot_points(3,:),'k'); %draw line
                end
            end
            hold off
        end
        function [pos_array, ind_array, z_array] = get_blade_pos(this)
            % 1) get the corners
            p{1} = this.HT * this.Blade.get_p1(); % in world frame
            p{2} = this.HT * this.Blade.get_p2();
            p{3} = this.HT * this.Blade.get_p3();
            p{4} = this.HT * this.Blade.get_p4();
            pos_temp = [p{1} p{2} p{3} p{4}];
            ind_array = get_coord(pos_temp,this.dx,this.minx,this.miny);
            for i = 1:4
                pi{i} =  ind_array(:,i);
            end
            
            % 2) get the lines            
            for l = 1:4 %from p1-p2, p2-p3, p3-p4, p4-p1
                if l<4 %just for cycling (p4-p1)
                    x1 = p{l}(1,1); y1 = p{l}(2,1);
                    x2 = p{l+1}(1,1); y2 = p{l+1}(2,1);
                elseif l==4
                    x1 = p{l}(1,1); y1 = p{l}(2,1);
                    x2 = p{1}(1,1); y2 = p{1}(2,1);
                end
                r = sqrt((x2-x1)^2 + (y2-y1)^2);
                inc = round((r-this.dx)/this.dx);
                p_ = zeros(3,1);
                for i = 1:inc
                    p_(1) = x1 + i*this.dx*(x2-x1)/r;
                    p_(2) = y1 + i*this.dx*(y2-y1)/r;
                    p_(3) = 1;
                    pi_ = get_coord(p_,this.dx,this.minx,this.miny);
                    ind_array = [ind_array pi_];
                end
            end
            
            % 3) get the area
            A = sortrows(ind_array')'; %sorted along first row
            % [1 1 1 2 2 3 3 3 
            %  1 2 3 1 3 1 2 3] ex small rectangle
            ind_temp = [];
            for i = A(1,1):A(1,end) % 1:3 one vertical line at a time
                temp = A(2,A(1,:)==i); % if we are looking at i = 2 temp = [1 3]
                for j = 1:temp(end)-temp(1)% j = 1:3-1 
                    if temp(1)+j < temp(end) && ~any(ismember(temp(1)+j,temp))
                        ind_temp = [ind_temp [i;temp(1)+j]]; %filling up
                    end
                end
            end
            ind_array = [ind_array ind_temp];
            pos_array = ind_array*this.dx + [this.minx-this.dx; this.miny-this.dx];
            pos_array = [pos_array;ones(1,size(pos_array,2))];
            pos_array_rel = this.get_invHT * pos_array;
            
            % 4) get corresponding z
            z_array = this.get_z(pos_array_rel);
%             z_array = -0.3*ones(1,size(pos_array,2)); %setting all depth
%             to be uniform
        end
        
        function z = get_z(this,p) %the depth of blade (it depends on where the position is)
            r = 0.5; %0.4;
            % this particular get_z is for
            % rectangle with x size of 1
            z = zeros(1,size(p,2));
            for i = 1:size(p,2)
                z(i) = -sqrt(r^2 - p(1,i)^2);
                if ~isreal(z(i))
                    z(i) = 0;
                end
            end
        end
        
        function [pos_array, ind_array] = get_displaced_pos(this)
            % 1) get the corners
            p{1} = this.HT * this.BladeFront.get_p1(); % in world frame
            p{2} = this.HT * this.BladeFront.get_p2();
            p{3} = this.HT * this.BladeFront.get_p3();
            p{4} = this.HT * this.BladeFront.get_p4();
            pos_temp = [p{1} p{2} p{3} p{4}];
            ind_array = get_coord(pos_temp,this.dx,this.minx,this.miny);
            for i = 1:4
                pi{i} =  ind_array(:,i);
            end
            % 2) get the lines            
            for l = 1:4 %from p1-p2, p2-p3, p3-p4, p4-p1
                if l<4 %just for cycling (p4-p1)
                    x1 = p{l}(1,1); y1 = p{l}(2,1);
                    x2 = p{l+1}(1,1); y2 = p{l+1}(2,1);
                elseif l==4
                    x1 = p{l}(1,1); y1 = p{l}(2,1);
                    x2 = p{1}(1,1); y2 = p{1}(2,1);
                end
                r = sqrt((x2-x1)^2 + (y2-y1)^2);
                inc = round((r-this.dx)/this.dx);
                p_ = zeros(3,1);
                for i = 1:inc
                    p_(1) = x1 + i*this.dx*(x2-x1)/r;
                    p_(2) = y1 + i*this.dx*(y2-y1)/r;
                    p_(3) = 1;
                    pi_ = get_coord(p_,this.dx,this.minx,this.miny);
                    ind_array = [ind_array pi_];
                end
            end
            % 3) get the area
            A = sortrows(ind_array')'; %sorted along first row
            % [1 1 1 2 2 3 3 3 
            %  1 2 3 1 3 1 2 3] ex small rectangle
            ind_temp = [];
            for i = A(1,1):A(1,end) % 1:3 one vertical line at a time
                temp = A(2,A(1,:)==i); % if we are looking at i = 2 temp = [1 3]
                for j = 1:temp(end)-temp(1)% j = 1:3-1 
                    if temp(1)+j < temp(end) && ~any(ismember(temp(1)+j,temp))
                        ind_temp = [ind_temp [i;temp(1)+j]]; %filling up
                    end
                end
            end
            ind_array = [ind_array ind_temp];
            pos_array = ind_array*this.dx + [this.minx-this.dx; this.miny-this.dx];          
        end
    end
%     methods (Static)
%         function plot_set(pos_cell)
%             for i = 1:length(pos_cell)
%                 pos_cell{i} =this.HT4 * this.PartsV{i}.pos_array;
%             end
%         end
%     end
end