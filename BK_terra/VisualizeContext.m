%%
classdef VisualizeContext < handle
    
    properties
        heightMatrix; %height matrix
%         hmapCellContext; %to keep a copy of HmapCellContext
        n; %grid per meter
        minx; miny;
        maxx; maxy; 
        resolution; %how many square boxes on one axis?
        yaw_angle_view = -60; 
        pitch_angle_view = 35;
        setupComplete = false;
        gridX; gridY; 
        hFigure; hAxes; hSurf; hLines = {}; 
    end
    properties (Constant)
        default_dT = 0.1;
        mycolormap = [
            0.5159    0.3224    0.2053  %my colormap for visualization
            0.5357    0.3348    0.2132
            0.5556    0.3472    0.2211
            0.5754    0.3596    0.2290
            0.5952    0.3720    0.2369
            0.6151    0.3844    0.2448
            0.6349    0.3968    0.2527
            0.6548    0.4092    0.2606
            0.6746    0.4216    0.2685
            0.6944    0.4340    0.2764
            0.7143    0.4464    0.2843
            0.7341    0.4588    0.2922
            0.7540    0.4712    0.3001
            0.7738    0.4836    0.3080
            0.7937    0.4960    0.3159
            0.8135    0.5084    0.3238
            0.8333    0.5208    0.3317
            0.8532    0.5332    0.3396
            0.8730    0.5456    0.3475
            0.8929    0.5580    0.3554
            0.9127    0.5704    0.3633
            0.9325    0.5828    0.3712
            0.9524    0.5952    0.3790
            0.9722    0.6076    0.3869
            0.9921    0.6200    0.3948
            1.0000    0.6324    0.4027
            1.0000    0.6448    0.4106
            1.0000    0.6572    0.4185
            1.0000    0.6696    0.4264
            1.0000    0.6820    0.4343
            1.0000    0.6944    0.4422
            1.0000    0.7068    0.4501
            1.0000    0.7192    0.4580
            1.0000    0.7316    0.4659
            1.0000    0.7440    0.4738
            1.0000    0.7564    0.4817
            1.0000    0.7688    0.4896
            1.0000    0.7812    0.4975];
    end
    %%
    methods
        % constructor synchronize with hmapCellContext at the beginning
        function this = VisualizeContext(minx,miny,maxx,maxy,n)
            this.n = n;
            this.minx = minx; this.miny = miny;
            this.maxx = maxx; this.maxy = maxy;
%             [hmap, ~] = hmapCellContext.get_Matrix;
%             this.heightMatrix = hmap;
%             this.resolution = hmapCellContext.resolution;
%             this.minx = hmapCellContext.minx; this.miny = hmapCellContext.miny;
%             this.maxx = hmapCellContext.maxx; this.maxy = hmapCellContext.maxy;
        end
        
        function visualize(this,hmatrix,lines)
            if nargin == 0
                hmatrix = this.heightMatrix; %if hmatrix is not given use what you have saved
            elseif nargin == 1
                this.heightMatrix = hmatrix;  %if no lines are given
            else
            end
            % Preprocedural Step for Drawing and Animating 
            if ~this.setupComplete
                this.setupComplete = true;
                % some preprocedure for visualization
                terrain_sizex = this.maxx - this.minx;
                xmin = this.minx;
                xmax = this.maxx;
                gridx = linspace(xmin,xmax,this.n*terrain_sizex+1)';
                gridX = [];
                % make gridX
                for i = 1:terrain_sizex*this.n+1
                    gridX = [gridX gridx];
                end
                
                terrain_sizey = this.maxy-this.miny;
                ymin = this.miny;
                ymax = this.maxy;
                gridy = linspace(ymin,ymax,this.n*terrain_sizey+1);
                gridY = [];
                % make gridY
                for i = 1:terrain_sizey*this.n+1
                    gridY = [gridY; gridy];
                end
                this.gridY = gridY;
                this.gridX = gridX;
                
                this.hFigure = figure;
                this.hAxes = axes(this.hFigure); 
                hold(this.hAxes, 'on');
                this.hSurf = surf(this.hAxes, this.gridX, this.gridY, hmatrix); %draw heightmap
                if nargin == 3 %lines given
                    for lineidx = 1:length(lines)
                        this.hLines{end+1} = plot3(this.hAxes, lines{lineidx}(1,:),lines{lineidx}(2,:),lines{lineidx}(3,:),'k-');
                    end
                end
                this.hSurf.EdgeColor = [0.5 0.5 0.5]; 
                colormap(this.mycolormap);
                axis(this.hAxes, [xmin xmax ymin ymax -xmax/5 xmax])
                xlabel(this.hAxes, 'x'); ylabel(this.hAxes, 'y'); zlabel(this.hAxes, 'z');
                this.hAxes.View = [this.yaw_angle_view this.pitch_angle_view];
                daspect([1 1 1]);
                hold(this.hAxes, 'off');
            else
                % Actual drawing happens here 
                [az,el] = view;
                hold(this.hAxes, 'on');
                this.hSurf.ZData = hmatrix; 

                if nargin == 3 %no lines given
                    for lineidx = 1:length(lines)
                        this.hLines{lineidx}.XData = lines{lineidx}(1,:); 
                        this.hLines{lineidx}.YData = lines{lineidx}(2,:);
                        this.hLines{lineidx}.ZData = lines{lineidx}(3,:);
%                         plot3(lines{lineidx}(1,:),lines{lineidx}(2,:),lines{lineidx}(3,:),'k-');
                    end
                end
                hold(this.hAxes, 'off');
                drawnow                
            end       
        end
        
        function SaveFigure(this, name)
            if nargin == 0
                name = 'figure'; 
            end 
            savefig(name); 
        end

    end
    
    methods(Static)
        function smethodsVisualize(hmapCellContext)            
            % some preprocedure for visualization
            terrain_sizex = hmapCellContext.maxx - hmapCellContext.minx;
            xmin = hmapCellContext.minx;
            xmax = hmapCellContext.maxx;
            gridx = linspace(xmin,xmax,hmapCellContext.resolution)';
            gridX = [];
            % make gridX
            for i = 1:hmapCellContext.resolution
                gridX = [gridX gridx];
            end
            
            terrain_sizey = hmapCellContext.maxy-hmapCellContext.miny;
            ymin = hmapCellContext.miny;
            ymax = hmapCellContext.maxy;
            gridy = linspace(ymin,ymax,hmapCellContext.resolution);
            gridY = []; 
            % make gridY
            for i = 1:hmapCellContext.resolution
                gridY = [gridY; gridy];
            end
            s = surf(gridX,gridY,hmapCellContext.get_Matrix); %draw heightmap
            s.EdgeColor = 'none';
            colormap(hmapCellContext.mycolormap);
            axis([xmin xmax ymin ymax -xmax/5 xmax])
            xlabel('x'); ylabel('y'); zlabel('z');
            drawnow
%             view([this.yaw_angle_view this.pitch_angle_view])
   
        end
    end
end