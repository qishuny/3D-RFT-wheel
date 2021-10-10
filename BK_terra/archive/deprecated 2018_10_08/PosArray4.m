%% give me pos array (particles) of an object
% it will have center pos x,y,z,th, and n particle points 
classdef PosArray4 < handle
    properties
        x; y; z; th;
        pos_array; %[x1 x2 ... xn;
                   % y1 y2 ... yn;
                   % z1 z2 ... zn
                   % 1  1  ... 1];
    end
    

    methods 
        
        function this = PosArray4(pos_array,x,y,z,th) %give pos array in 3xn 
            this.x = x; %where its center is located at
            this.y = y;
            this.z = z;
            this.th = th;
            this.pos_array = get_HT4(x,y,z,th)...
                *[pos_array;ones(1,size(pos_array,2))];
            
        end
    end
    
end