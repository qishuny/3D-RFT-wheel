%% index object for converting particle points to indices and making sure none of them overlap
classdef IndexHandler < handle
    
    properties

        minx; miny;
        maxx; maxy; 
        resolution; %how many square boxes on one axis?
        dx; %each increment is how many meters?
    end

    %%
    methods
        % constructor synchronize with hmapCellContext at the beginning
        function this = IndexHandler(dx,minx,miny,maxx,maxy)
            this.dx = dx;
            this.minx = minx; this.miny = miny;
            this.maxx = maxx; this.maxy = maxy;
            
%             this.resolution = hmapCellContext.resolution;
%             this.dx = hmapCellContext.get_dx();
%             this.minx = hmapCellContext.minx; this.miny = hmapCellContext.miny;
%             this.maxx = hmapCellContext.maxx; this.maxy = hmapCellContext.maxy;
        end
        
        function [indices_norepeat,z_vec] = getIndex(this,pointcloud) %get indices with no repeat
            x_vec = pointcloud(1,:) - this.minx;
            y_vec = pointcloud(2,:) - this.miny;
            x_idx = round(x_vec/this.dx)+1;
            y_idx = round(y_vec/this.dx)+1;
            indices = [x_idx; y_idx];
            [indices_norepeat,IA,~] = unique(indices','rows');
            indices_norepeat = indices_norepeat';
            z_vec = pointcloud(3,IA);
            
        end
             

    end
end