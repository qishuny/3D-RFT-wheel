%%
classdef Part < handle
    properties
        % contour: for now restrict to rectangles
        minx; miny;
        maxx; maxy;
        % f(x,y) equation for depth 
    end
    

    methods 
        
        function this = Part(minx,miny,maxx,maxy,x,y)
            % the midpoint is at x,y
            % size of the rectangle is determined by
            % minx miny maxx maxy
            if nargin == 0
                this.minx = -0.3; this.miny = -0.3;
                this.maxx = 0.3;  this.maxy = 0.3;
            elseif nargin == 6
                this.minx = minx + x; this.miny = miny + y;
                this.maxx = maxx + x; this.maxy = maxy + y;
            else
                error("number of inputs must be appropriate")
            end
        end
        function set_corners(this,p1,p3)
            this.minx = p1(1);
            this.miny = p1(2);
            this.maxx = p3(1);
            this.maxy = p3(2);
        end
        function p1 = get_p1(this)
            p1 = [this.minx; this.miny; 1];
        end
        function p2 = get_p2(this)
            p2 = [this.maxx; this.miny; 1];
        end
        function p3 = get_p3(this)
            p3 = [this.maxx; this.maxy; 1];
        end
        function p4 = get_p4(this)
            p4 = [this.minx; this.maxy; 1];
        end
    end
    
end