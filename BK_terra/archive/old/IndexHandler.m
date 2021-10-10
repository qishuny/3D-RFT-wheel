%%
classdef IndexHandler < handle %handle class to avoid unnecessary copying in function calls
    %use copy() if necessary
    
    properties (GetAccess = 'public', SetAccess = 'private')
        dx;
        minx;
        miny; %from Hmap properties
    end
    
    methods
        function this = IndexHandler(Hmap) %IndexHandler and Hmap must be synchronous
            if ~isa(Hmap,'Hmap')
                error('Input must be class:Hmap');
            end
            this.dx = Hmap.get_dx;
            this.minx = Hmap.minx;
            this.miny = Hmap.miny;
        end
        
        function [pos_array, ind_array, z_array] = get_contact(Vehicle)
            if ~isa(Vehicle,'Vehicle');
                error('Input must be class:Vehicle');
            end
            % 1) get the corners
            p{1} = Vehicle.HT * Vehicle.Blade.get_p1(); % in world frame
            p{2} = Vehicle.HT * Vehicle.Blade.get_p2();
            p{3} = Vehicle.HT * Vehicle.Blade.get_p3();
            p{4} = Vehicle.HT * Vehicle.Blade.get_p4();
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
        
        function displaced_ind = get_displaced(Vehicle,vel)
            if ~isa(Vehicle,'Vehicle')
                error('1st Input must be class:Vehicle');
            end
            [pos_array, ~,~] = get_contact(Vehicle);
            vel_unit = vel/(vel(1)^2+vel(2)^2);
            displaced_pos = zeros(size(pos_array));
            for i = 1:length(displaced_pos)
                displaced_pos(1,i) = pos_array(1,i) + vel_unit(1)*this.dx;
                displaced_pos(2,i) = pos_array(2,i) + vel_unit(2)*this.dx;
            end
            displaced_ind = get_coord(displaced_pos,this.dx,this.minx,this.miny);
        end
        
        
        
    end
end