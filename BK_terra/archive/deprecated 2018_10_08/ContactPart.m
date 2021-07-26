%%
classdef ContactPart < handle
    properties
        x; y; th; %relative to Vehicle
        lines;
        
        dx;
        points;
        points_idx;
    end
    

    methods 
        
        function this = ContactPart(x,y,th,dx,lines)
            if nargin == 0
                this.x = 0;
                this.y = 0;
                this.th = 0;
                this.dx = 0;
                % right angle triangle
                this.lines{1} = Line(0,0,1,0); %x1 y1 x2 y2
                this.lines{2} = Line(1,0,0,1);
                this.lines{3} = Line(1,0,0,0);
            elseif nargin == 4
                this.x = x;
                this.y = y;
                this.th = th;
                this.dx = dx;
                this.lines{1} = Line(0,0,1,0); %x1 y1 x2 y2
                this.lines{2} = Line(1,0,0,1);
                this.lines{3} = Line(1,0,0,0);
            elseif nargin == 5
                this.x = x;
                this.y = y;
                this.th = th;
                this.dx = dx;
                if ~isa(lines,'cell')
                    error('Edges input must be cell array');
                end
                this.lines = lines; %must be cell array
            else
                error("Number of inputs must be four or five (x,y,th,dx,lines)")
            end
        end
        
        function HT = getHT(this)
            HT = [cos(this.th) -sin(this.th) this.x;
                sin(this.th) cos(this.th) this.y
                0 0 1];
        end
        function invHT = getinvHT(this)
            R = [cos(this.th) -sin(this.th);
                sin(this.th) cos(this.th)];
            p = -R'*[this.x;this.y];
            invHT = [cos(this.th) sin(this.th) p(1)
                -sin(this.th) cos(this.th) p(2)
                0 0 1];
        end
        function relpos =  get_pos(this)
            relpos = [];
        end
        
        function area_idx = get_coord_inside(this)       
            if ~isa(this.lines,'cell')
                error('Input must be cell array of Line s');
            end
            dx = this.lines{1}.dx;
            edges_idx = [];
            for i = 1:3
                edges_idx = [edges_idx get_coord(this.lines{i}.points,dx)];
            end
            edges_idx_sorted = sortrows(edges_idx')';
            area_idx = edges_idx; %initialize with line indices
            for x_idx = min(edges_idx_sorted(1,:)):max(edges_idx_sorted(1,:))
                % from left to right
                temp = edges_idx_sorted(:,edges_idx_sorted(1,:)==x_idx);
                % points corresponding to x_idx == 1, .. 2,3.. so on
                % [2 2      or     [1 1 1 ... 1
                %  1 10]            1 2 3 ... 11]
                y_idx = min(temp(2,:));
                while y_idx < max(temp(2,:))
                    % points between (x_idx==1,y_min) and (x_idx==1,y_max)
                    if ~ismember(y_idx,temp(2,:))
                        area_idx = horzcat(area_idx,[x_idx;y_idx]);
                    end
                    y_idx = y_idx + 1;
                end
            end
        end
        
        
    end
    
end