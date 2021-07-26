%% -- initialize --
clear
close all
clc

% for Desktop
% cd('C:\Users\BKIM\Documents\GitHub\terramech\BK_terra');
% for XPS 
cd('C:\Users\BrianKim\Documents\GitHub\terramech\BK_terra');

disp('Directory Changed');


%% -- ToDO: 1) Fix blocking flow 2) Update only part that matters 3) Get optimized dt parameter

% parameter for wheel 
% diameter: 0.096m
% width of wheel: 0.048m
% angle of repose 29 degrees

%% -- 3d plot practice --
close all

% circle with respect to center of object at (0,0)
% translate all this with HT
r = 0.5;
t = 0:pi/10:2*pi;
st = sin(t);
ct = cos(t);
uno = ones(size(t));
x_points = ct;
width = 0.4;
y_points = width*uno;
z_points = st;
plot_points{1} = [x_points;y_points;z_points;ones(size(x_points))];
y_points = -width*uno;
plot_points{2} = [x_points;y_points;z_points;ones(size(x_points))];
for i = 1:length(t)
    plot_points{2+i} = [ct(i), ct(i);-width,width;st(i),st(i)];
end
plot_points{end+1} = [0 0;0 2;0.7 0.7];
plot_points{end+1} = [3 3;0 2;0.7 0.7];
plot_points{end+1} = [0 3;0 0;0.7 0.7];
plot_points{end+1} = [0 3;2 2;0.7 0.7];


%draw 
figure
hold on
for i = 1:length(plot_points)
    plot3(plot_points{i}(1,:),plot_points{i}(2,:),plot_points{i}(3,:),'k');
end
view([-50 30]);
% plot3([0 0],[0 2],[0.7 0.7],'black');
% plot3([3 3],[0 2],[0.7 0.7],'black');



% figure
% hold on
% % plot3(ct,0*z,st,'black')
% plot3(ct,-0.4*uno,st,'black')
% plot3(ct,0.4*uno,st,'black')
% % plot3(0.8*ct,-0.32*uno,0.8*st,'k--')
% % plot3(0.8*ct,0.4*uno,0.8*st,'k--')
% for i = 1:length(t)
%     plot3([ct(i) ct(i)],[-0.4 0.4],[st(i) st(i)],'k');
% end
% xlabel('x'); ylabel('y'); zlabel('z');
% view([-50 30]);
% % hold off
% 
% % rectangle
% x1 = 0:0.1:3;
% y1 = 0*ones(size(x1));
% z1 = 0.7*ones(size(x1));
% uno = ones(size(x1));
% plot3(x1,y1,z1,'black');
% plot3(x1,2*uno,z1,'black');
% plot3([0 0],[0 2],[0.7 0.7],'black');
% plot3([3 3],[0 2],[0.7 0.7],'black');

%% % ------------ 1D testing ------------------
close all
clear
clc

% using matrix.. the form I am comfortable with
% 1 Impulse
impulse_area = 1;

height = zeros(1,100);
alpha = 1;
dt = 0.005; %0.004;
dx = 0.1;
impulse_height = impulse_area/dx;
height(1,end) = impulse_height;

% height = zeros(1,10);
% alpha = 1;
% dt = 0.4;
% dx = 1;
% impulse_height = impulse_area/dx;
% height(1,end) = impulse_height;

% dx = 0.5;
% maxx = 10;
% resolution = round(maxx/dx);
% height = zeros(1,resolution);
% alpha = 1;
% dt = 0.1;
% impulse_height = impulse_area/dx;
% height(1,end) = impulse_height;

% dx = 0.01;
% maxx = 10;
% resolution = round(maxx/dx);
% height = zeros(1,resolution);
% alpha = 1;
% dt = 0.00005;
% impulse_height = impulse_area/dx;
% height(1,end) = impulse_height;

angle_repose = 0.5;

xaxis = zeros(size(height));
for i = 1:length(xaxis)
    xaxis(i) = i*dx;
end

figure
plot(height)
axis([1 10 0 10])
temp_data = [];
pause
for i = 1:3000 %100  %300
    height = simple_find_new_height(height,dx,alpha,dt,angle_repose);
    pause(0.02)
    plot(xaxis,height)
    axis([1 10 0 10])
    if mod(i,10) == 0
        temp_data = [temp_data; height];
    end
end

data{10} = temp_data;


%% 2 Medium apart impulse
for impulse_location = 1:9
    height = zeros(1,10);
    height(1,10) = 5;
    height(1,impulse_location) = 5;
    figure
    plot(height)
    temp_data = [];
    for i = 1:100
        alpha = 0.1;
        dt = 1;
        dx = 1;
        height = simple_find_new_height(height,dx,alpha,dt);
        pause(0.02)
        plot(height)
        axis([1 10 0 5])
        if mod(i,10) == 0
            temp_data = [temp_data; height];
        end
    end
    data{impulse_location} = temp_data;
end

%% plotting different steady states at 
figure
subplot(5,2,1);
plot(data{1}(end,:));
for i = 2:10
    subplot(5,2,i);
    plot(data{i}(end,:));
end
title('Steady state of different impulse responses');

%% Impulse after one goes steady state
for next_impulse = 1:10
    height = zeros(1,10);
    alpha = 0.1;
    dt = 1;
    dx = 1;
    height(1,10) = 5;
    figure
    plot(height)
    temp_data = [];
    for i = 1:100
        height = simple_find_new_height(height,dx,alpha,dt);
        pause(0.01)
        plot(height)
        axis([1 10 0 5])
        if mod(i,10) == 0
            temp_data = [temp_data; height];
        end
    end
    % after ss add right next
    height(1,next_impulse) = height(1,next_impulse) + 5;
    for i = 1:100
        height = simple_find_new_height(height,dx,alpha,dt);
        pause(0.01)
        plot(height)
        axis([1 10 0 5])
        if mod(i,10) == 0
            temp_data = [temp_data; height];
        end
    end
    data_afterss{next_impulse} = temp_data;
end
%% ------------- Cell Method Practice ------------------
% practicing cell method 
% use cell method?? yea for practice
clear
close all
clc

% initialize 
left_corner_x = 0;
left_corner_y = 0;
right_corner_x = 10;
right_corner_y = 10;
resolution = 101;
hmapCellContext = HmapCellContext(left_corner_x,left_corner_y,right_corner_x,right_corner_y,resolution);
visualizeContext = VisualizeContext(hmapCellContext); %synchronize them first
indexHandler = IndexHandler(hmapCellContext); %synchronize with IndexHandler
% check initial state
hmap = hmapCellContext.get_Matrix();
visualizeContext.visualize(hmap);

% perturb the system 
perturb_location = [5;5];
indices = indexHandler.getIndex(perturb_location)
hmapCellContext.HmapCell_array{indices(1),indices(2)}.height = 10;
dT = 0.1;
doVis = true;
hmapCellContext.steadyState(dT,doVis);
[hmap,~] = hmapCellContext.get_Matrix();

visualizeContext.visualize(hmap);

% s = surf(gridX,gridY,hmatrix);

% % expand(beginning_cell);
% begin_cell = HmapCell();
% begin_cell.i = 1; begin_cell.j = 1;
% begin_cell.x = 0; begin_cell.y = 0;
% 
% %setup
% grid_array = cell(10,10);
% for i = 1:10
%     for j = 1:10
%         grid_array{i,j} = HmapCell(i,j,i*0.1-0.1,j*0.1-0.1);
%         grid_array{i,j}.height = 0;
%     end
% end
% 
% %set back pointers for neighbors
% for i = 1:10
%     for j = 1:10
%         if i>1 
%             grid_array{i,j}.neighbors{end+1} = grid_array{i-1,j};
%         end
%         if j>1
%             grid_array{i,j}.neighbors{end+1} = grid_array{i  ,j-1};
%         end
%         if j<10
%             grid_array{i,j}.neighbors{end+1} = grid_array{i  ,j+1};
%         end
%         if i<10
%             grid_array{i,j}.neighbors{end+1} = grid_array{i+1,j};
%         end
%     end
% end
% 
% % % % update based on angle of repose
% % % q_sand = sign(slope)*(abs(slope)-thres)*(abs(slope)>thres)*0.001/dx*dt;
% % % 
% % % update = sign(slope).*(abs(slope) - thres).*(abs(slope)>thres)*adx;
% % % alphadx = 0.001 / dx * dt;
% % % 
% % % h_grid = grid_array{4,3};
% % % slope = (h_grid.neighbors{1}.height - h_grid.height)/dx;
% % % q_sand = k*dx*sign(slope)*(abs(slope)-angle_repose)*(abs(slope)>angle_repose);
% % % area = dx*dx;
% % % delta_h = q_sand*dt/area;
% % % h_temp = h + delta_h;
% 
% % update to an impulse response
% dx = 0.1;
% dt = 0.1;
% k = 0.01;
% angle_repose = 0.6;
% 
% grid_array{5,5}.height = 10;
% 
% for i = 1:size(grid_array,1)
%     for j = 1:size(grid_array,2)
%         sand_from_neighbors = 0;
%         h_grid = grid_array{i,j}; %copy of it
%         for neigh_idx = 1:length(grid_array{i,j}.neighbors)
%             neighbor = h_grid.neighbors{neigh_idx};
%             slope = (neighbor.height - h_grid.height)/dx;
%             sand_from_neighbors = sand_from_neighbors + ...
%                 k*dx*sign(slope)*(abs(slope)-angle_repose)*(abs(slope)>angle_repose);
%         end
%         delta_h = sand_from_neighbors*dt/(dx*dx);
%         % save new height on height_temp, update all at once later
%         grid_array{i,j}.height_temp = grid_array{i,j}.height + delta_h;
%     end
% end
% 
% % update all at once
% for i = 1:size(grid_array,1)
%     for j = 1:size(grid_array,2)
%         grid_array{i,j}.height = grid_array{i,j}.height_temp;
%     end
% end
% 
% h_map = zeros(size(grid_array));
% for i = 1:size(grid_array,1)
%     for j = 1:size(grid_array,2)
%         h_map(i,j) = grid_array{i,j}.height;
%     end
% end
% h_map


%% -- 2.5D representation testing ---
pos = Bull.get_blade_pos;
sortrows(pos')';
figure
plot(pos(1,:),pos(2,:),'.');
axis equal
size(pos)

Bull.move(0,1);
pos = Bull.get_blade_pos;
figure
plot(pos(1,:),pos(2,:),'.');
axis equal
size(pos)

Bull.move(0,2);
pos = Bull.get_blade_pos;
figure
plot(pos(1,:),pos(2,:),'.');
axis equal
size(pos)


%% -- Hmap Sand Rolling Demo --
clear 
close all
n = 10; %grid per m
minx = -8;
miny = -8;
maxx = 8;
maxy = 8;
Sand = Hmap(minx,miny,maxx,maxy,n);
Sand.dT = 2;

Sand.visualize
pause
Sand.add_sand(-8,8,-20);
for i = 1:40
Sand.add_sand(i*Sand.get_dx,0,10)
end
Sand.visualize
pause
Sand.steady_state('on');
for j = 1:30
    Sand.add_sand(4,j*Sand.get_dx,10);
end
Sand.visualize
pause
Sand.steady_state('on');
for i = 1:40
    Sand.add_sand(-4+i*Sand.get_dx,-4,-5);
end
Sand.visualize
pause
Sand.steady_state('on');


%% -- PoseContext practice ---
close all
clear
Body = PoseContext([],0,0,0,0);
line = [0 0;0 0;-1 1;1 1];
Wheel1_axis = PoseContext([],1,1,0,0);
Wheel1_axis.parent = Body;
Body.child = Wheel1_axis;
Wheel1 = PoseContext(line,0,0,0,0);
Wheel1.parent = Wheel1_axis;
Wheel1_axis.child = Wheel1;

p_world = Body.HT*Wheel1_axis.HT*Wheel1.HT*Wheel1.particles;

figure
% plot3(line(1,:),line(2,:),line(3,:),'k');
plot3(p_world(1,:),p_world(2,:),p_world(3,:),'k');
xlabel('x'); ylabel('y'); zlabel('z');

hold on
for i = 1:8
    pause(0.1)
Wheel1.rotate_x(pi/8);
p_world = Body.HT*Wheel1_axis.HT*Wheel1.HT*Wheel1.particles;
plot3(p_world(1,:),p_world(2,:),p_world(3,:),'k');
end
hold off
% Wheel1.rotate_x(pi/8);
% p_world = Body.HT*Wheel1_axis.HT*Wheel1.HT*Wheel1.particles;
% plot3(p_world(1,:),p_world(2,:),p_world(3,:),'k');


%% ---- Impuse Sand Rolling For Experiment V2 --
clear 
close all

gridsizes = [5 10 20 40];
heights = cell(1,length(gridsizes));
for gridsizeidx = 1:length(gridsizes)
    n = gridsizes(gridsizeidx); %grid per m
    minx = -1;
    miny = -1;
    maxx = 1;
    maxy = 1;
    Sand = Hmap(minx,miny,maxx,maxy,n);
    % Sand.dT = 2;
    
    %-----------initialize visualizeContext
    visualizeContext = VisualizeContext(Sand.minx,Sand.miny,Sand.maxx,Sand.maxy,Sand.n);
    
    % ----------initialize indexHandler
    indexHandler = IndexHandler(Sand.get_dx,Sand.minx,Sand.miny,Sand.maxx,Sand.maxy);
    
    %-----------visualize hmap
    visualizeContext.visualize(Sand.matrix);
    
    % ----------perturb
    impulse_volume = 0.1;
    impulse_height = impulse_volume/Sand.get_dx^2;
    Sand.add_sand(0,0,impulse_height)
    
    % ---------update
    Sand.update_onestep;
    pointsidx  = 1;

    dataper = 100;
    initial = 10;
    maximumheight = zeros(1,2*round(1/Sand.get_dx^2));
    for i = 1:2*round(1/Sand.get_dx^2)
        Sand.update_onestep;

        if mod(i,dataper) == 0
            fromx = -1;
            tox = 1;
            points_from = indexHandler.getIndex([fromx;0;0]);
            points_to = indexHandler.getIndex([tox;0;0]);
            
            sliceofpoints{pointsidx} = Sand.matrix(points_from(1):points_to(1),points_from(2):points_to(2))';
            xpoints{pointsidx} = fromx:Sand.get_dx:tox;
            heightmaps{pointsidx} = Sand.matrix;
            pointsidx = pointsidx + 1;
        end
        maximumheight(1,i) = max(max(Sand.matrix));
    end
%     visualizeContext.visualize(Sand.matrix);
    steps = 1:i;
    maximumheight = [maximumheight;steps];
    heights{gridsizeidx} = maximumheight;
    
%     % line with tand(29) slope
%     angleofrepose = 29;
%     x0 = 0;
%     y0 = maximumheight(end);
%     x_repose = linspace(-y0/tand(angleofrepose) + x0, x0);
%     height_repose = tand(angleofrepose)*(x_repose-x0) + y0;
%     
%     
%     figure
%     hold on
%     for i = 1:length(sliceofpoints)
%         plot(xpoints{i},sliceofpoints{i},'Color',[0.6151    0.3844    0.2448],'DisplayName',sprintf('%d',dataper*(i)))
%     end
%     plot(x_repose,height_repose,'b','MarkerSize',9,'DisplayName','Angle of Repose')
%     hold off
%     axis equal
  
end
figure
hold on
for i = 1:length(heights)
loglog(heights{i}(2,:),heights{i}(1,:))
end
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off
legend('n=5','n=10','n=20','n=40')
%% ---- Hmap Sand Rolling Demo V2 --
clear 
close all
n = 10; %grid per m
minx = -8;
miny = -8;
maxx = 8;
maxy = 8;
Sand = Hmap(minx,miny,maxx,maxy,n);
% Sand.dT = 2;

%-----------initialize visualizeContext
visualizeContext = VisualizeContext(Sand.minx,Sand.miny,Sand.maxx,Sand.maxy,Sand.n);

% ----------initialize indexHandler
indexHandler = IndexHandler(Sand.get_dx,Sand.minx,Sand.miny,Sand.maxx,Sand.maxy);

%-----------visualize hmap
% visualizeContext.visualize(Sand.matrix);

tic
% ----------perturb and visualize
Sand.add_sand(-8,8,-20); %perturb
for i = 1:40
Sand.add_sand(i*Sand.get_dx,0,10)
end
visualizeContext.visualize(Sand.matrix);
pause

Sand.update_onestep;
while ~Sand.reachedSteadystate
    Sand.update_onestep;
    visualizeContext.visualize(Sand.matrix);
end


%perturb
for j = 1:30
    Sand.add_sand(4,j*Sand.get_dx,10);
end
visualizeContext.visualize(Sand.matrix);
pause
%update and visualize
Sand.update_onestep;
while ~Sand.reachedSteadystate
    Sand.update_onestep;
    visualizeContext.visualize(Sand.matrix);
end


%perturb 
for i = 1:40
    Sand.add_sand(-4+i*Sand.get_dx,-4,-5);
end
visualizeContext.visualize(Sand.matrix);
pause
%update and visualize
Sand.update_onestep;
while ~Sand.reachedSteadystate
    Sand.update_onestep;
    visualizeContext.visualize(Sand.matrix);
end
toc


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

angle = [0 22.5 45 66.5 90];
depth = [0.011, 0.012];
experimentresults = cell(1,length(angle)*length(depth));
expidx = 1;

for depthidx = 2 %1:length(depth)
    for angleidx = 4:5 %1:length(angle)
        n = 200; %grid per m; make it at most 0.02
        minx = -0.4;
        miny = -0.4;
        maxx = 0.4;
        maxy = 0.4;
        Sand = Hmap(minx,miny,maxx,maxy,n);
        
        
        %----------------setting up body position
        bodyxpos = 0;
        bodyypos = 0 - 0.1;
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
        final_pos = [bodyxpos bodyypos+0.2 bodyzpos]';
        final_vel = [0 initial_vel(2) 0]';
        final_acc = [0 0 0]';
        tfinal = norm(final_pos - initial_pos)/norm(initial_vel); %10;
        dT = 1;
        % gives me vector of X,Y,Z for trajectory
        [X,Y,Z,dX,dY,dZ,theta] = trajectoryContext.minJerk(initial_pos,initial_vel,initial_acc,final_pos,final_vel,final_acc,tfinal,dT);
        
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
        saveas(gcf,sprintf('wheelexperiment_depth%dmm_angle%0.1fdeg.png',1000*depth(depthidx),angle(angleidx)));
        
        experimentresults{expidx} = data;
        expidx = expidx + 1;
    end
end
        

%% -- Hmap Wheel Demo---
clear 
close all
n = 10; %grid per m
minx = -8;
miny = -8;
maxx = 8;
maxy = 8;
Sand = Hmap(minx,miny,maxx,maxy,n);
Sand.dT = 2;

Bull = Vehicle(0,0,0,Sand.get_dx);

% set HSmap for my wheel
% 1) the x,y position of each points
r = 0.5;  %radius of the wheel
xgrid = -r:Sand.get_dx/2:r; %particles' xpos
ygrid = -0.2:Sand.get_dx/2:0.2; %2; %particles' ypos
xy_vec = []; 
for i = 1:length(ygrid)
    temp = [xgrid; ones(size(xgrid))*ygrid(i)];
    xy_vec = [xy_vec temp]; % [x1 x2 x3 ... xn;
                                  %  y1 y2 y3 ... yn];
end
% pos_array = [pos_array; ones(1,size(pos_array,2))];
% 2) the z position of each points
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
%wheel 1
wheel1 = PosArray3(xy_vec,z_vec,1.5,2.5,pi/2);
% PartVisualization(wheel1_PosArray3,'wheel1');
Bull.PosArray3_cell{1} = wheel1;

wheel_x = 1.5; wheel_y = -2.5; wheel_th = -pi/4; %wheel 2
Bull.PosArray3_cell{2} = PosArray3(xy_vec,z_vec,wheel_x,wheel_y,wheel_th);
wheel_x = -1.5; wheel_y = -2.5; wheel_th = pi/2; %wheel 3
Bull.PosArray3_cell{3} = PosArray3(xy_vec,z_vec,wheel_x,wheel_y,wheel_th);
wheel_x = -1.5; wheel_y = 2.5; wheel_th = pi/2; %wheel 4
Bull.PosArray3_cell{4} = PosArray3(xy_vec,z_vec,wheel_x,wheel_y,wheel_th);

%for visualization : using PoseContext
% Body = PoseContext([],Bull.get_x(),Bull.get_y(),0,Bull.get_th());
% line = [0 0;0 0;-1 1;1 1];
% Wheel1_axis = PoseContext([],1,1,0,0);
% Wheel1_axis.parent = Body;
% Body.child = Wheel1_axis;
% Wheel1 = PoseContext(line,0,0,0,0);
% Wheel1.parent = Wheel1_axis;
% Wheel1_axis.child = Wheel1;
% 
% p_world = Body.HT*Wheel1_axis.HT*Wheel1.HT*Wheel1.particles;


% getting the indices of the wheel 
blade_idx_all = []; z_array_all = [];
for i = 1:length(Bull.PosArray3_cell) %for four wheels
    PosArray3_pos_world = Bull.HT*Bull.PosArray3_cell{i}.pos_array;
%     blade_idx = get_coord(Bull.PosArray3_cell{i}.pos_array,Sand.get_dx,Sand.minx,Sand.miny);
    blade_idx = get_coord(PosArray3_pos_world,Sand.get_dx,Sand.minx,Sand.miny);
    [blade_idx_norepeat,IA,~] = unique(blade_idx','rows');
    blade_idx_norepeat = blade_idx_norepeat';
    z_array_norepeat = z_vec(IA);
    blade_idx_all = [blade_idx_all blade_idx_norepeat];
    z_array_all = [z_array_all z_array_norepeat];
end

Bull.synchronize(Sand);
Bull.move(0,0);
Sand.steady_state('on');

% circle with respect to center of object at (0,0)
% translate all this with HT
t = 0:pi/10:2*pi;
st = r*sin(t);
ct = r*cos(t);
uno = ones(size(t));
x_points = ct;
width = 0.2;
y_points = width*uno;
z_points = st;
plot_points{1} = [x_points;y_points;z_points]; %circle line 1
y_points = -width*uno;
plot_points{2} = [x_points;y_points;z_points]; %circle line 2
for i = 1:length(t)
    plot_points{2+i} = [ct(i), ct(i);-width,width;st(i),st(i)]; %tire lines
end

for i = 1:length(plot_points) %wheel 1
    posarray_cell{i} = PosArray4(plot_points{i},1.5,2.5,0,pi/2);
end
Bull.PartsV{1} = PartVisualization(posarray_cell,'wheel1');
posarray_cell = {};
for i = 1:length(plot_points) % wheel 2
    posarray_cell{i} = PosArray4(plot_points{i},1.5,-2.5,0,-pi/4);
end
Bull.PartsV{2} = PartVisualization(posarray_cell,'wheel2');
posarray_cell = {};
for i = 1:length(plot_points) % wheel 3
    posarray_cell{i} = PosArray4(plot_points{i},-1.5,-2.5,0,pi/2);
end
Bull.PartsV{3} = PartVisualization(posarray_cell,'wheel3');
posarray_cell = {};
for i = 1:length(plot_points) % wheel 4
    posarray_cell{i} = PosArray4(plot_points{i},-1.5,2.5,0,pi/2);
end
Bull.PartsV{4} = PartVisualization(posarray_cell,'wheel4');


% for idx_wheels = 1:4
%     Bull.PartsV{i} = 
%     for i = 1:length(plot_points)
%         Bull.PartsV{i} = ...
%             PosArray4(plot_points{i},1.5,2.5,0,pi/2);
% %         PosArray4_set = PosArray4(wheel_plot{i},1.5,2.5,0,pi/2); 
%     end
% end

% plot_cell = Bull.get_plot_pos();
% hold on
% for i = 1:length(plot_cell)
%     plot3(plot_cell{i}(1,:),plot_cell{i}(2,:),plot_cell{i}(3,:),'k');
% end
% hold off
% 
% pause
% Sand.resolve_collision(blade_idx_norepeat,z_array_norepeat,[1;0]);
% Sand.steady_state('on');

% use the data from the indices to update the Sandmap accordingly
Sand.resolve_collision(blade_idx_all,z_array_all,[0;1]);
Sand.steady_state('off');



pause
hmap_overtime = cell(1,25);
for i = 1:25
    Bull.move(1,0);
    blade_idx_all = []; z_array_all = [];
    % get the blade indices of all wheels
    for j = 1:4
        PosArray3_pos_world = Bull.HT*Bull.PosArray3_cell{j}.pos_array;
        %     blade_idx = get_coord(Bull.PosArray3_cell{i}.pos_array,Sand.get_dx,Sand.minx,Sand.miny);
        blade_idx = get_coord(PosArray3_pos_world,Sand.get_dx,Sand.minx,Sand.miny);
        [blade_idx_norepeat,IA,~] = unique(blade_idx','rows');
        blade_idx_norepeat = blade_idx_norepeat';
        z_array_norepeat = z_vec(IA);
        blade_idx_all = [blade_idx_all blade_idx_norepeat];
        z_array_all = [z_array_all z_array_norepeat];
    end
    Sand.resolve_collision(blade_idx_all,z_array_all,[0;1]);
    Sand.steady_state('off');
    hold on    
    Bull.visualize();
    hold off
    drawnow;
%     plot_cell = Bull.get_plot_pos();
%     hold on
%     for k = 1:length(plot_cell)
%         plot3(plot_cell{k}(1,:),plot_cell{k}(2,:),plot_cell{k}(3,:),'k');
%     end

%     pause
    
%     hold off
%     pause

    hmap_overtime{i} = Sand.matrix; %for later visualization
end

% later visualization (faster)
% figure
% pause
% for i = 1:length(hmap_overtime)
%     pause(0.2);
%     Hmap.visualize_matrix(hmap_overtime{i},Sand.minx,Sand.miny,Sand.maxx,Sand.maxy,Sand.n);
% end

%% -- Hmap Testing ---
clear
% close all

n = 10; %grid per m
minx = -8;
miny = -8;
maxx = 8;
maxy = 8;
Sand = Hmap(minx,miny,maxx,maxy,n);
Sand.dT = 2;
% Blade = Part(-0.4,0.2,  0.4,0.4);
% Blade = Part(0.2,-0.4,  0.6,0.4);
Blade = Part(-0.5,-0.2,  0.5,0.2, 0.5,0.5); %minx miny maxx maxy x y
Blade_f = Part(-0.5,-0.2,  0.5,0.2, 0.5,0.5+0.4);

% Blade2 = Part(-1.5,0.2,  -0.5,0.6);
% Blade2_f = Part(-1.5,0.8,  -0.5,1);

Bull = Vehicle(0,0,0,Sand.get_dx,Blade,Blade_f);


% set HSmap for my wheel
% 1) the position of each points
r = 0.5;
xgrid = -r:Sand.get_dx/2:r;
ygrid = -0.2:Sand.get_dx/2:0.2;
xy_vec = [];
for i = 1:length(ygrid)
    temp = [xgrid; ones(size(xgrid))*ygrid(i)];
    xy_vec = [xy_vec temp];
end
% pos_array = [pos_array; ones(1,size(pos_array,2))];
% 2) the z depth of each points
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
wheel_x = 1;
wheel_y = 1;
wheel1V = PosArray3(xy_vec,z_vec,wheel_x,wheel_y);
Vehicle.PosArray3 = wheel1V;

blade_idx = get_coord(Vehicle.PosArray3.pos_array,Sand.get_dx,Sand.minx,Sand.miny);
[blade_idx_norepeat,IA,~] = unique(blade_idx','rows');
blade_idx_norepeat = blade_idx_norepeat';
z_array_norepeat = z_vec(IA);

% set PosArray3 for wheel front
% 1) set position
r = 0.5;
xgrid = -r:Sand.get_dx/2:r;
ygrid = -0.2:Sand.get_dx/2:0.2;
xy_vec = [];
for i = 1:length(ygrid)
    temp = [xgrid; ones(size(xgrid))*ygrid(i)];
    xy_vec = [xy_vec temp];
end
wheel_x = 1;
wheel_y = 1.4;
wheel1_displaced_PosArray3 = PosArray3(xy_vec,z_vec,wheel_x,wheel_y);
Vehicle.PosArray3_displaced = wheel1_displaced_PosArray3;

displaced_idx = get_coord(Vehicle.PosArray3_displaced.pos_array,Sand.get_dx,Sand.minx,Sand.miny);
[displaced_idx_norepeat,IA,~] = unique(displaced_idx','rows');
displaced_idx_norepeat = displaced_idx_norepeat';






Bull.synchronize(Sand);
Bull.move(0,0);

Sand.steady_state('on');
% Sand.add_sand(0,0,10);
% Sand.steady_state('on');

% Sand.sync(blade_idx_norepeat,displaced_idx_norepeat,z_array_norepeat);
% Sand.steady_state('on');

Sand.resolve_collision(blade_idx_norepeat,z_array_norepeat,[0;1]);
Sand.steady_state('on');

% [pos, b_ind, z] = Bull.get_blade_pos();
% [d_pos, d_ind] = Bull.get_displaced_pos();
% Sand.sync(b_ind,d_ind,z);
% Sand.steady_state('on');

% Sand.add_sand(1,1,10);
% Sand.add_sand(-1,-1,10);
% % Sand.add_sand(0,0.8,10);
% Sand.steady_state('on');

% for i = 1:6
%     Bull.move(0,-1);
%     [pos, b_ind, z] = Bull.get_blade_pos();
%     [d_pos, d_ind] = Bull.get_displaced_pos();
%     Sand.sync(b_ind,d_ind,z);
%     Sand.steady_state('off');
% end

% for i = 1:15
%     Bull.move(1,0);
%     [pos, b_ind, z] = Bull.get_blade_pos();
%     [d_pos, d_ind] = Bull.get_displaced_pos();
%     Sand.sync(b_ind,d_ind,z);
%     Sand.steady_state('off');
% end
