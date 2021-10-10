%%
classdef Line < handle
    properties
        x1;
        y1;
        x2;
        y2;
        points;
        dx;
    end
    

    methods 
        
        function this = Line(x1,y1,x2,y2,dx)
            if nargin == 0
                this.x1 = 0;
                this.y1 = 0;
                this.x2 = 1;
                this.y2 = 1;
                this.dx = 0.1;
                set_dx(this,this.dx);
            elseif nargin == 4
                this.x1 = x1;
                this.y1 = y1;
                this.x2 = x2;
                this.y2 = y2;
                this.dx = 0.1;
                set_dx(this,this.dx);
            elseif nargin == 5
                this.x1 = x1;
                this.y1 = y1;
                this.x2 = x2;
                this.y2 = y2;
                this.dx = dx;
                set_dx(this,this.dx);
            else 
                error("Number of inputs must be five (x1,y1,x2,y2,dx)")
            end
        end
        
        function set_begin(this,x,y)
            this.x1 = x;
            this.y1 = y;
        end
        
        function set_end(this,x,y) 
            this.x2 = x;
            this.y2 = y;
        end
       
        function set_dx(this,dx)

            y1 = this.y1; x2 = this.x2; %#ok<*PROPLC>
            x1  = this.x1; y2 = this.y2;
            r = sqrt((x2-x1)^2+(y2-y1)^2);
            inc = round(r/dx);
            for i = 1:inc
                this.points(1,i) = x1 + i*dx*(x2-x1)/r;
                this.points(2,i) = y1 + i*dx*(y2-y1)/r;
                this.points(3,i) = 1;
            end

        end
    end
   
end