%%
classdef HSmap < handle
    properties
        x; y; th;
        pos_array; %[x1 x2 ... xn;
                   % y1 y2 ... yn;
                   % 1  1  ... 1];
        z_array;
    end
    

    methods 
        
        function this = HSmap(pos_array,z_array,x,y,th) %pos array and z array with respect to 0,0 
            this.x = x; %where its center is located at
            this.y = y;
            this.th = th;
            this.pos_array = get_HT(x,y,th)*[pos_array;ones(1,size(pos_array,2))];
            this.z_array = z_array;
            
        end
    end
    
end