classdef HmapCell < handle
    
    properties
        height = 0;       
        height_temp = 0;
%         left_cell;
%         right_cell;
        
        neighbors;
        x; y; i; j;
        
        isActive; %flag for telling us if this needs updating

    end
    %%
    methods
        function this = HmapCell(i,j,x,y)
            if nargin == 0
                %             this.height = h;
                this.neighbors = {};
                %             this.idx = i;
            elseif nargin == 4
                this.i = i; this.j = j;
                this.x = x; this.y = y;
            end
                
        end
        
        function flag = isSame(this,element) %for use in Queue
            % two HmapCell is the same if their i,j are the same
            flag = false;
            if this.i == element.i && this.j == element.j
                flag = true;
            end
        end

    end
end
