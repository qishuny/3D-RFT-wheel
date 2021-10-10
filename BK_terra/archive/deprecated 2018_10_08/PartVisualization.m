%% object for combining PosArray, name for visualization
% wheel1, wheel2, wheel3, ... mainly for easier organization
classdef PartVisualization < handle
    properties
        lines; %will contain cell of particles(PosArray4)
        name;
    end
    

    methods 
        
        function this = PartVisualization(posarraycell,str) %pos array and z array with respect to 0,0 
            this.name = str;
            this.lines = posarraycell;
        end
        
        function cell = get_posArray(this)
            cell = this.lines;
        end
    end
    
end