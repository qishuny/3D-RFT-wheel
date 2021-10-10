%% -- wheel experiment for comparison---
clear 
close all

% parameter for wheel 
% diameter: 0.096m
% width of wheel: 0.048m
% angle of repose 29 degrees

% params for different experiments 
% angle = 0, 22.5, 45, 66.5, 90
% depth = 0.011, 0.012

angle = [0 22.5 45 66.5 90]; % degrees 
depth = [0.01, 0.02, 0.03]; % meters
experimentresults = cell(1,length(angle)*length(depth));
expidx = 1;
DISTANCE = 0.6; % how far wheel move in each settings (don't go over maxx - minx) 


for depthidx = 1:length(depth)
    for angleidx = 1:length(angle)
        n = 200; %grid per m; make it at most 0.02
        minx = -0.4;
        miny = -0.4;
        maxx = 0.4;
        maxy = 0.4;
        Sand = Hmap(minx,miny,maxx,maxy,n);
        
        
        %----------------setting up body position
        bodyxpos = 0;
        bodyypos = 0 - DISTANCE/2;
        bodyzpos = 0.096/2-depth(depthidx);
        bodytheta = 0; %body is straigh, change wheel angle
        Body = Vehicle(bodyxpos,bodyypos,bodytheta,Sand.get_dx);
        BodyPose = PoseContext([],bodyxpos,bodyypos,bodyzpos,bodytheta);
        
        
        %----------------setting pointcloud for wheel
        % set particlepoints for my wheel points wrt origin of wheel pose
        % 1) the x,y position of each points
        x_width = 0.096/2;  %1/2 width of blade
        xgrid = -x_width:Sand.get_dx/2:x_width; %particles' xpos
        y_length = 0.048/2; %1/2 length of blade just enough to cover like 2~3 grids points
        ygrid = -y_length:Sand.get_dx/2:y_length; %particles' ypos
        xy_vec = [];
        for i = 1:length(ygrid)
            temp = [xgrid; ones(size(xgrid))*ygrid(i)];
            xy_vec = [xy_vec temp]; % [x1 x2 x3 ... xn;
            %  y1 y2 y3 ... yn];
        end
        % 2) the z position of each points
        r = x_width;
        z_vec = zeros(1,length(xy_vec));
        for i = 1:length(z_vec)
            temp = -sqrt(r^2 - xy_vec(1,i)^2);
            if isreal(temp)
                z_vec(i) = temp; %circle with r
            else
                z_vec(i) = 0;
                disp('z is imaginary');
            end
        end
        particlepoints = [xy_vec;z_vec;ones(size(z_vec))];
        
        % join two poses: body and wheel
        wheelxpos = 0; %always relative position of pose to parent
        wheelypos = 0;
        wheelzpos = 0;
        wheeltheta = (90 - angle(angleidx))*pi/180;
        WheelPose = PoseContext(particlepoints,wheelxpos,wheelypos,wheelzpos,wheeltheta);
        WheelPose.parent = BodyPose;
        BodyPose.child = WheelPose;
        
        
        
        % --------------- setting lines for visualization
        % wheel, we need many lines to draw wheel
        t = 0:pi/10:2*pi;
        st = r*sin(t);
        ct = r*cos(t);
        uno = ones(size(t));
        
        %circle line 1
        x_points = ct;
        y_points = y_length*uno;
        z_points = st;
        lines{1} = [x_points;y_points;z_points;uno]; %circle line 1
        %circle line 2
        y_points = -y_length*uno;
        lines{2} = [x_points;y_points;z_points;uno]; %circle line 2
        %tire lines
        for i = 1:length(t)
            lines{2+i} = [ct(i), ct(i);-y_length,y_length;st(i),st(i);1 1]; %tire lines
        end
        %WheelPose now contains particlespoints for collision and
        %lines for visualization
        WheelPose.lines = lines;
        
        %-----------initialize visualizeContext
        visualizeContext = VisualizeContext(Sand.minx,Sand.miny,Sand.maxx,Sand.maxy,Sand.n);
        
        % ----------initialize indexHandler
        indexHandler = IndexHandler(Sand.get_dx,Sand.minx,Sand.miny,Sand.maxx,Sand.maxy);
        
        % ----------initialize TrajectoryContext
        trajectoryContext = TrajectoryContext();
        
        % ----------making trajectory
        initial_pos = [bodyxpos bodyypos bodyzpos]';
        initial_vel = [0 Sand.get_dx/2 0]';
        initial_acc = [0 0 0]';
        final_pos = [bodyxpos bodyypos+DISTANCE bodyzpos]';
        final_vel = [0 initial_vel(2) 0]';
        final_acc = [0 0 0]';
        tfinal = norm(final_pos - initial_pos)/norm(initial_vel); %10;
        dT = 1;
        % gives me vector of X,Y,Z for trajectory
        [X,Y,Z,dX,dY,dZ,theta] = trajectoryContext.minJerk(initial_pos,initial_vel,initial_acc,...
                                                            final_pos,final_vel,final_acc,tfinal,dT);
        
        %-----------Resolve collision, visualize blade, the initial step
        % Bull.synchronize(Sand);
        
        particlepointsWorld = BodyPose.HT4*WheelPose.HT4*WheelPose.particles;
        [indices, zvec] = indexHandler.getIndex(particlepointsWorld);
        velocity = [0;1];
        Sand.blocked_idx = [];
        Sand.resolve_collision(indices,zvec,velocity);
        
        lines_wcoord = cell(size(lines)); %getting the world coordinate lines for visualization
        for i = 1:length(WheelPose.lines)
            lines_wcoord{i} = BodyPose.HT4*WheelPose.HT4*WheelPose.lines{i};
        end
        
        Sand.update_onestep;
        for i=1:140
            Sand.update_onestep;
            %     visualizeContext.visualize(Sand.matrix,lines_wcoord);
        end
%         visualizeContext.visualize(Sand.matrix,lines_wcoord);

        % --------------- THE MAIN SIMULATION ------------------------
        %------------Resolve collision, visualize blade, following steps
        savedataSand = cell(1,length(X));
        savedatalines = cell(1,length(X));
        for i = 1:length(X)
            % body's trajectory
            r = sqrt(dX(i)^2+dY(i)^2);
            BodyPose.HT4 = [dY(i)/r dX(i)/r 0 X(i);
                -dX(i)/r dY(i)/r 0 Y(i);
                0 0 1 Z(i);
                0 0 0 1];
            
            particlepointsWorld = BodyPose.HT4*WheelPose.HT4*WheelPose.particles;
            [indices, zvec] = indexHandler.getIndex(particlepointsWorld);
            velocity = [dX(i),dY(i)];
            
            %     Sand.Steadystate(indices,zvec,velocity);
            Sand.blocked_idx = [];
            Sand.resolve_collision(indices,zvec,velocity);
            %     Sand.steady_state('off');
            Sand.update_onestep;
            Sand.resolve_collision(indices,zvec,velocity);
            for steadystate=1:140
                Sand.update_onestep;
                Sand.resolve_collision(indices,zvec,velocity);
                %         Sand.resolve_collision(indices,zvec,velocity);
                %         visualizeContext.visualize(Sand.matrix,lines_wcoord);
            end
            
            lines_wcoord = cell(size(lines)); %getting the world coordinate lines for visualization
            for j = 1:length(WheelPose.lines)
                lines_wcoord{j} = BodyPose.HT4*WheelPose.HT4*WheelPose.lines{j};
            end
%             visualizeContext.visualize(Sand.matrix,lines_wcoord);
            %     savedataSand{i} = Sand.matrix;
            %     savedatalines{i} = lines_wcoord;
        end
        
        %visualizing the last step
        lines_wcoord = cell(size(lines)); %getting the world coordinate lines for visualization
        for j = 1:length(WheelPose.lines)
            lines_wcoord{j} = BodyPose.HT4*WheelPose.HT4*WheelPose.lines{j};
        end
        visualizeContext.visualize(Sand.matrix,lines_wcoord);
        
        % which points am I going to sample for comparison?
        points_from = indexHandler.getIndex([-0.2;0;0]);
        points_to = indexHandler.getIndex([0.2;0;0]);
        
        sliceofpoints = Sand.matrix(points_from(1):points_to(1),points_from(2):points_to(2));
        sliceofpoints = sliceofpoints';
        xpoints = -0.2:Sand.get_dx:0.2;
        data.points = [xpoints;sliceofpoints];
        data.wheelangle = angle(angleidx);
        data.wheeldepth = depth(depthidx);
        data.sand = Sand;
        data.lines = lines_wcoord;
%         saveas(gcf,sprintf('wheelexperiment_depth%dmm_angle%0.1fdeg.png',1000*depth(depthidx),angle(angleidx)));
        visualizeContext.SaveFigure(sprintf('wheelexp_depth%d_angle%0.1fdeg.fig',1000*depth(depthidx),angle(angleidx)));
        experimentresults{expidx} = data;
        expidx = expidx + 1;
    end
end

       