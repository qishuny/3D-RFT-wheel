%% general pose object for handling series of poses
classdef PoseContext < handle
    properties
        HT; HT4; parent = []; child = [];
        particles; %[x1 x2 ... xn;   for collision check
                   % y1 y2 ... yn;
                   % z1 z2 ... zn
                   % 1  1  ... 1];
        lines = {}; %cell of lines for visualization
    end
    

    methods 
        
        function this = PoseContext(particlePoints,x,y,z,th) %x,y,z,th with respect to its parent
            this.HT4 = get_HT4(x,y,z,th);
            this.HT = get_HT(x,y,th);
            this.particles= particlePoints;
        end
        
        function rotate_z(this,th)
            this.HT = get_HT(0,0,th)*this.HT;
            this.HT4 = get_HT4(0,0,0,th)*this.HT4;
        end
        function x = get_x(this)
            x = this.HT4(1,4);
        end
        function y = get_y(this)
            y = this.HT4(2,4);
        end
        function z = get_z(this)
            z = this.HT4(3,4);
        end
        function th = get_th(this)
            th = atan2(this.HT4(2,1),this.HT4(2,2)); 
        end
        
    end
    
end