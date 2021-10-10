%% -- One Wheel Demo---
clear 
close all

% parameter for wheel 
% diameter: 0.096m
% width of wheel: 0.048m
% angle of repose 29 degrees

n = 100; %grid per m; make it at most 0.02
minx = -0.4;
miny = -0.4;
maxx = 0.4;
maxy = 0.4;
Sand = Hmap(minx,miny,maxx,maxy,n);


%----------------setting up body position
bodyxpos = 0;
bodyypos = miny + abs(miny/5);
bodyzpos = 0;
bodytheta = 0;
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
wheeltheta = pi/2;
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
final_pos = [bodyxpos bodyypos+(maxy-initial_pos(2))-maxy/3 bodyzpos]';
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
Sand.resolve_collision(indices,zvec,velocity);

Sand.steady_state('off'); %we will no longer use visualization from Hmap class
lines_wcoord = cell(size(lines)); %getting the world coordinate lines for visualization
for i = 1:length(WheelPose.lines)
    lines_wcoord{i} = BodyPose.HT4*WheelPose.HT4*WheelPose.lines{i};
end

Sand.update_onestep;
for i=1:140
    Sand.update_onestep;
%     visualizeContext.visualize(Sand.matrix,lines_wcoord);
end
visualizeContext.visualize(Sand.matrix,lines_wcoord);
pause
%------------Resolve collision, visualize blade, following steps
savedataSand = cell(1,length(X));
savedatalines = cell(1,length(X));
for i = 1:length(X)
    % body's trajectory
    r = sqrt(dX(i)^2+dY(i)^2);
    BodyPose.HT4 = [dY(i)/r dX(i)/r 0 X(i);
        -dX(i)/r dY(i)/r 0 Y(i);
        0 0 1 0;
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
    for steadystate=1:200
        Sand.update_onestep;
        Sand.resolve_collision(indices,zvec,velocity);
%         Sand.resolve_collision(indices,zvec,velocity);
%         visualizeContext.visualize(Sand.matrix,lines_wcoord);
    end

    lines_wcoord = cell(size(lines)); %getting the world coordinate lines for visualization
    for j = 1:length(WheelPose.lines)
        lines_wcoord{j} = BodyPose.HT4*WheelPose.HT4*WheelPose.lines{j};
    end
    visualizeContext.visualize(Sand.matrix,lines_wcoord);
    savedataSand{i} = Sand.matrix;
    savedatalines{i} = lines_wcoord;
end

points_from = indexHandler.getIndex([-0.2;0;0]);
points_to = indexHandler.getIndex([0.2;0;0]);

sliceofpoints = Sand.matrix(points_from(1):points_to(1),points_from(2):points_to(2));
xpoints = -0.2:Sand.get_dx:0.2;
figure
plot(xpoints,sliceofpoints)
axis equal

% figure
% pause
% for i = 1:length(savedataSand)
%     visualizeContext.visualize(savedataSand{i},savedatalines{i});
% end


% 
% %% -- Hmap Wheel Demo---
% clear 
% close all
% n = 10; %grid per m
% minx = -8;
% miny = -8;
% maxx = 8;
% maxy = 8;
% Sand = Hmap(minx,miny,maxx,maxy,n);
% Sand.dT = 2;
% 
% Bull = Vehicle(0,0,0,Sand.get_dx);
% 
% % set HSmap for my wheel
% % 1) the position of each points
% r = 2; %0.5;  %radius of the wheel
% xgrid = -r:Sand.get_dx/2:r; %particles' xpos
% ygrid = -0.5:Sand.get_dx/2:0.5; %2; %particles' ypos
% pos_array = []; 
% for i = 1:length(ygrid)
%     temp = [xgrid; ones(size(xgrid))*ygrid(i)];
%     pos_array = [pos_array temp]; % [x1 x2 x3 ... xn;
%                                   %  y1 y2 y3 ... yn];
% end
% % pos_array = [pos_array; ones(1,size(pos_array,2))];
% % 2) the z depth of each points
% z_array = zeros(1,length(pos_array));
% for i = 1:length(z_array)
%     temp = -sqrt(r^2 - pos_array(1,i)^2);
% =======
% %%
% clear
% close all
% clc
% 
% % initialize 
% left_corner_x = 0;
% left_corner_y = 0;
% right_corner_x = 10;
% right_corner_y = 10;
% resolution = 101;
% 
% hmapCellContext = HmapCellContext(left_corner_x,left_corner_y,right_corner_x,right_corner_y,resolution);
% 
% x_idx = round(size(hmapCellContext.HmapCells,1)/2);
% y_idx = round(size(hmapCellContext.HmapCells,2)/2);
% impulse_vol = 1;
% impulse_height = impulse_vol/hmapCellContext.get_dx()^2;
% hmapCellContext.HmapCells{x_idx,y_idx}.height = impulse_height;
% % hmapCellContext.HmapCells{5,5}.isActive = true;
% 
% doVis = true;
% tic
% hmapCellContext.steadyState(0.15,doVis)
% toc
% 
% 
% % set pointcloud for my wheel
% % 1) the x,y,z position of each points
% radius = 0.5;  %radius of the wheel
% xpoints = -radius:hmapCellContext.get_dx/2:radius; %particles' xpos
% % ypoints = -0.2:hmapCellContext.get_dx/2:0.2; %particles' ypos
% ypoints = zeros(size(xpoints));
% zpoints = zeros(size(xpoints));
% for i = 1:length(zpoints)
%     temp = -sqrt(radius^2 - xpoints(i)^2);
% >>>>>>> 4bfdfdea83c3832cf1b9c37b19b8d4d8e67e6414
%     if isreal(temp)
%         zpoints(i) = temp; %circle with r
%     else
%         zpoints(i) = 0;
%         disp('z is imaginary');
%     end
% end
% <<<<<<< HEAD
% %wheel 1
% % wheel1 = PosArray3(pos_array,z_array,2.5,2.5,pi/2);
% % PartVisualization(wheel1_PosArray3,'wheel1');
% % Bull.PosArray3_cell{1} = wheel1;
% 
% wheel_x = 4.5; wheel_y = -2.5; wheel_th = -pi/4; %wheel 2
% Bull.PosArray3_cell{end+1} = PosArray3(pos_array,z_array,wheel_x,wheel_y,wheel_th);
% wheel_x = -4.5; wheel_y = -2.5; wheel_th = pi/2; %wheel 3
% Bull.PosArray3_cell{end+1} = PosArray3(pos_array,z_array,wheel_x,wheel_y,wheel_th);
% % wheel_x = -1.5; wheel_y = 2.5; wheel_th = pi/2; %wheel 4
% % Bull.PosArray3_cell{4} = PosArray3(pos_array,z_array,wheel_x,wheel_y,wheel_th);
% 
% %for visualization : using PoseContext
% Body = PoseContext([],Bull.get_x(),Bull.get_y(),0,Bull.get_th());
% line = [0 0;0 0;-1 1;1 1];
% Wheel1_axis = PoseContext([],1,1,0,0);
% =======
% pointcloud = [xpoints;ypoints;zpoints];
% bull_xpos = 2;
% bull_ypos = 2;
% bull_th = 0;
% 
% Body = PoseContext([],bull_xpos,bull_ypos,bull_th);
% Wheel1_axis = PoseContext([],1,1,0,0); %relative position to body
% >>>>>>> 4bfdfdea83c3832cf1b9c37b19b8d4d8e67e6414
% Wheel1_axis.parent = Body;
% Body.child = Wheel1_axis;
% Wheel1 = PoseContext(pointcloud,0,0,0,0); %relative position to wheelaxis
% Wheel1.parent = Wheel1_axis;
% Wheel1_axis.child = Wheel1;
% 
% <<<<<<< HEAD
% p_world = Body.HT*Wheel1_axis.HT*Wheel1.HT*Wheel1.particles;
% 
% 
% % getting the indices of the wheel 
% blade_idx_all = []; z_array_all = [];
% for i = 1:length(Bull.PosArray3_cell) %for four wheels
%     PosArray3_pos_world = Bull.HT*Bull.PosArray3_cell{i}.pos_array;
% %     blade_idx = get_coord(Bull.PosArray3_cell{i}.pos_array,Sand.get_dx,Sand.minx,Sand.miny);
%     blade_idx = get_coord(PosArray3_pos_world,Sand.get_dx,Sand.minx,Sand.miny);
%     [blade_idx_norepeat,IA,~] = unique(blade_idx','rows');
%     blade_idx_norepeat = blade_idx_norepeat';
%     z_array_norepeat = z_array(IA);
%     blade_idx_all = [blade_idx_all blade_idx_norepeat];
%     z_array_all = [z_array_all z_array_norepeat];
% end
% 
% Bull.synchronize(Sand);
% Bull.move(0,0);
% Sand.steady_state('on');
% 
% % circle with respect to center of object at (0,0)
% % translate all this with HT
% t = 0:pi/10:2*pi;
% st = r*sin(t);
% ct = r*cos(t);
% uno = ones(size(t));
% x_points = ct;
% width = 0.2;
% y_points = width*uno;
% z_points = st;
% plot_points{1} = [x_points;y_points;z_points];
% y_points = -width*uno;
% plot_points{2} = [x_points;y_points;z_points];
% for i = 1:length(t)
%     plot_points{2+i} = [ct(i), ct(i);-width,width;st(i),st(i)];
% end
% 
% for i = 1:length(plot_points) %wheel 1
%     posarray_cell{i} = PosArray4(plot_points{i},1.5,2.5,0,pi/2);
%     %         PosArray4_set = PosArray4(wheel_plot{i},1.5,2.5,0,pi/2);
% end
% Bull.PartsV{1} = PartVisualization(posarray_cell,'wheel1');
% posarray_cell = {};
% for i = 1:length(plot_points) % wheel 2
%     posarray_cell{i} = PosArray4(plot_points{i},1.5,-2.5,0,0);
% end
% Bull.PartsV{2} = PartVisualization(posarray_cell,'wheel2');
% posarray_cell = {};
% for i = 1:length(plot_points) % wheel 3
%     posarray_cell{i} = PosArray4(plot_points{i},-1.5,-2.5,0,pi/2);
% end
% Bull.PartsV{3} = PartVisualization(posarray_cell,'wheel3');
% posarray_cell = {};
% for i = 1:length(plot_points) % wheel 4
%     posarray_cell{i} = PosArray4(plot_points{i},-1.5,2.5,0,pi/2);
% end
% Bull.PartsV{4} = PartVisualization(posarray_cell,'wheel4');
% 
% 
% % for idx_wheels = 1:4
% %     Bull.PartsV{i} = 
% %     for i = 1:length(plot_points)
% %         Bull.PartsV{i} = ...
% %             PosArray4(plot_points{i},1.5,2.5,0,pi/2);
% % %         PosArray4_set = PosArray4(wheel_plot{i},1.5,2.5,0,pi/2); 
% %     end
% % end
% 
% % plot_cell = Bull.get_plot_pos();
% % hold on
% % for i = 1:length(plot_cell)
% %     plot3(plot_cell{i}(1,:),plot_cell{i}(2,:),plot_cell{i}(3,:),'k');
% % end
% % hold off
% % 
% % pause
% % Sand.resolve_collision(blade_idx_norepeat,z_array_norepeat,[1;0]);
% % Sand.steady_state('on');
% 
% % use the data from the indices to update the Sandmap accordingly
% Sand.resolve_collision(blade_idx_all,z_array_all,[0;1]);
% Sand.steady_state('off');
% 
% 
% 
% pause
% hmap_overtime = cell(1,25);
% for i = 1:45
%     Bull.move(1,0);
%     blade_idx_all = []; z_array_all = [];
%     % get the blade indices of all wheels
%     for j = 1:length(Bull.PosArray3_cell)
%         PosArray3_pos_world = Bull.HT*Bull.PosArray3_cell{j}.pos_array;
%         %     blade_idx = get_coord(Bull.PosArray3_cell{i}.pos_array,Sand.get_dx,Sand.minx,Sand.miny);
%         blade_idx = get_coord(PosArray3_pos_world,Sand.get_dx,Sand.minx,Sand.miny);
%         [blade_idx_norepeat,IA,~] = unique(blade_idx','rows');
%         blade_idx_norepeat = blade_idx_norepeat';
%         z_array_norepeat = z_array(IA);
%         blade_idx_all = [blade_idx_all blade_idx_norepeat];
%         z_array_all = [z_array_all z_array_norepeat];
%     end
%     Sand.resolve_collision(blade_idx_all,z_array_all,[0;1]);
%     Sand.steady_state('off');
%     
% %     Bull.visualize();
%     
% %     plot_cell = Bull.get_plot_pos();
% %     hold on
% %     for k = 1:length(plot_cell)
% %         plot3(plot_cell{k}(1,:),plot_cell{k}(2,:),plot_cell{k}(3,:),'k');
% %     end
%     
% %     hold off
% %     pause
%     
%     hmap_overtime{i} = Sand.matrix; %for later visualization
% end
% 
% 
% figure
% pause
% for i = 1:length(hmap_overtime)
%     pause(0.07);
%     Hmap.visualize_matrix(hmap_overtime{i},Sand.minx,Sand.miny,Sand.maxx,Sand.maxy,Sand.n);
% end

% %%
% clear 
% close all
% n = 10; %grid per m
% minx = -8;
% miny = -8;
% maxx = 8;
% maxy = 8;
% Sand = Hmap(minx,miny,maxx,maxy,n);
% Sand.dT = 2;
% 
% Bull = Vehicle(0,0,0,Sand.get_dx);
% 
% % set HSmap for my wheel
% % 1) the position of each points
% r = 0.5;  %radius of the wheel
% xgrid = -r:Sand.get_dx/2:r; %particles' xpos
% ygrid = -0.2:Sand.get_dx/2:0.2; %particles' ypos
% pos_array = []; 
% for i = 1:length(ygrid)
%     temp = [xgrid; ones(size(xgrid))*ygrid(i)];
%     pos_array = [pos_array temp]; % [x1 x2 x3 ... xn;
%                                   %  y1 y2 y3 ... yn];
% end
% % pos_array = [pos_array; ones(1,size(pos_array,2))];
% % 2) the z depth of each points
% z_array = zeros(1,length(pos_array));
% for i = 1:length(z_array)
%     temp = -sqrt(r^2 - pos_array(1,i)^2);
%     if isreal(temp)
%         z_array(i) = temp; %circle with r
%     else
%         z_array(i) = 0;
%         disp('z is imaginary');
%     end
% end
% %wheel 1
% wheel1 = PosArray3(pos_array,z_array,1.5,2.5,pi/2);
% % PartVisualization(wheel1_PosArray3,'wheel1');
% Bull.PosArray3_cell{1} = wheel1;
% 
% wheel_x = 1.5; wheel_y = -2.5; wheel_th = -pi/4; %wheel 2
% Bull.PosArray3_cell{2} = PosArray3(pos_array,z_array,wheel_x,wheel_y,wheel_th);
% wheel_x = -1.5; wheel_y = -2.5; wheel_th = pi/2; %wheel 3
% Bull.PosArray3_cell{3} = PosArray3(pos_array,z_array,wheel_x,wheel_y,wheel_th);
% wheel_x = -1.5; wheel_y = 2.5; wheel_th = pi/2; %wheel 4
% Bull.PosArray3_cell{4} = PosArray3(pos_array,z_array,wheel_x,wheel_y,wheel_th);
% 
% %for visualization : using PoseContext
% Body = PoseContext([],Bull.get_x(),Bull.get_y(),Bull.get_z(),Bull.get_th());
% line = [0 0;0 0;-1 1;1 1];
% Wheel1_axis = PoseContext([],1,1,0,0);
% Wheel1_axis.parent = Body;
% Body.child = Wheel1_axis;
% Wheel1 = PoseContext(line,0,0,0,0);
% Wheel1.parent = Wheel1_axis;
% Wheel1_axis.child = Wheel1;
% 
% p_world = Body.HT*Wheel1_axis.HT*Wheel1.HT*Wheel1.particles;
% % r = 0.4;
% % dx = 0.1;
% % min_x = -r+dx/4;
% % max_x = r-dx/4;
% % x = min_x:dx/2:max_x;
% % z = zeros(size(x));
% % for i=1:length(x)
% %     z(i) = -sqrt(r^2 - x(i)^2);
% % end
% % plot(x,z)
% =======
% pointcloud_world = Body.HT*Wheel1_axis.HT*Wheel1.HT*Wheel1.particles;
% >>>>>>> 4bfdfdea83c3832cf1b9c37b19b8d4d8e67e6414


%% Checking how handle works
% clear
% close all
% clc
% 
% % initialize 
% left_corner_x = 0;
% left_corner_y = 0;
% right_corner_x = 10;
% right_corner_y = 10;
% resolution = 11;
% 
% hmapCellContext = HmapCellContext(left_corner_x,left_corner_y,right_corner_x,right_corner_y,resolution);
% 
% hmapCellContext.HmapCells{5,5}.height = 10;
% hmapCellContext.HmapCells{5,5}.isActive = true;
% 
% doVis = true;
% hmapCellContext.steadyState(0.1,doVis)
% hmapCellContext.activeQueue
% 




% hmapCellContext.active_HmapCells = {};
% hmapCellContext.active_HmapCells{end+1} = hmapCellContext.HmapCellarray{5,5};
% for neighidx = 1:length(hmapCellContext.HmapCellarray{5,5}.neighbors)
%     hmapCellContext.active_HmapCells{end+1} = hmapCellContext.HmapCellarray{5,5}.neighbors{neighidx};
%     hmapCellContext.HmapCellarray{5,5}.neighbors{neighidx}.isActive = true;
% end

%for activeCells, updateOnestep();
% temp_activeCells = hmapCellContext.active_HmapCells;






% dx = this.get_dx();
% k = this.flow_rate;
% angle_repose = this.angle_repose;
% 
% 
% % update only where the flow is happening (where previous slope
% % was over the threshold)
% for actidx = 1:length(this.activeHmapCells)
%     h_cell = this.activeHmapCells{actidx}; %copy of it
%     sand_from_neighbors = 0; %this will contain flow from all neighbors for this h_grid
%     
%     % sand flow q from all neighbors
%     for neigh_idx = 1:length(h_cell.neighbors)
%         neighbor = h_cell.neighbors{neigh_idx};
%         slope = (neighbor.height - h_cell.height)/dx;
%         if slope > angle_repose
%             neighbor.active = true;
%         end
%         sand_from_neighbors = sand_from_neighbors + ...
%             k*dx*sign(slope)*(abs(slope)-angle_repose)*(abs(slope)>angle_repose);
%     end
%     % sand flow rate q * dt = sand flow volume
%     delta_h = sand_from_neighbors*dt/(dx*dx);
%     % save new height on height_temp, update all at once later
%     this.activeHmapCells{i}.height_temp = this.activeHmapCells{i}.height + delta_h;
% end
% 
% % update all at once
% [hmap, hmap_temp] = this.get_Matrix();
% this.maxheightchange_laststep = max(max(hmap_temp - hmap)); %save how much changed for flag
% % update all at once
% for actidx = 1:length(this.activeHmapCells)
%     this.activeHmapCells{actidx}.height = this.activeHmapCells{actidx}.height_temp;
% end

% 
% 
% % initialize 
% left_corner_x = 0;
% left_corner_y = 0;
% right_corner_x = 10;
% right_corner_y = 10;
% resolution = 11;
% hmapCellContext = HmapCellContext(left_corner_x,left_corner_y,right_corner_x,right_corner_y,resolution);
% visualizeContext = VisualizeContext(hmapCellContext); %synchronize them first
% % check initial state
% hmap = hmapCellContext.get_Matrix();
% visualizeContext.visualize(hmap);
% 
% % perturb the system 
% hmapCellContext.HmapCell_array{5,5}.height = 10;
% hmapCellContext.steadyState();
% [hmap,~] = hmapCellContext.get_Matrix();
% 
% visualizeContext.visualize(hmap);

% clear
% close all
% clc
% 
% % expand(beginning_cell);
% begin_cell = Cell();
% begin_cell.i = 1; begin_cell.j = 1;
% begin_cell.x = 0; begin_cell.y = 0;
% 
% %setup
% grid_array = cell(10,10);
% for i = 1:10
%     for j = 1:10
%         grid_array{i,j} = Cell(i,j,i*0.1-0.1,j*0.1-0.1);
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
% % update to an impulse response
% dx = 0.1;
% dt = 0.1;
% k = 1;
% angle_repose = 0.6;
% 
% grid_array{5,5}.height = 10;
% 
% for i = 1:size(grid_array,1)
%     for j = 1:size(grid_array,2)
%         sand_from_neighbors = 0;
%         h_grid = grid_array{i,j}; %copy of it
%         for neigh_idx = 1:length(grid_array{3,3}.neighbors)
%             neighbor = h_grid.neighbors{neigh_idx};
%             slope = (neighbor.height - h_grid.height)/dx;
%             sand_from_neighbors = sand_from_neighbors + ...
%                 k*dx*sign(slope)*(abs(slope)-angle_repose)*(abs(slope)>angle_repose);
%         end
%         delta_h = sand_from_neighbors*dt/(dx*dx);
%         grid_array{i,j}.height_temp = grid_array{i,j}.height + delta_h;
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
% 
% 
% % %%
% % clear 
% % close all
% % n = 10; %grid per m
% % minx = -8;
% % miny = -8;
% % maxx = 8;
% % maxy = 8;
% % Sand = Hmap(minx,miny,maxx,maxy,n);
% % Sand.dT = 2;
% % 
% % Bull = Vehicle(0,0,0,Sand.get_dx);
% % 
% % % set HSmap for my wheel
% % % 1) the position of each points
% % r = 0.5;  %radius of the wheel
% % xgrid = -r:Sand.get_dx/2:r; %particles' xpos
% % ygrid = -0.2:Sand.get_dx/2:0.2; %particles' ypos
% % pos_array = []; 
% % for i = 1:length(ygrid)
% %     temp = [xgrid; ones(size(xgrid))*ygrid(i)];
% %     pos_array = [pos_array temp]; % [x1 x2 x3 ... xn;
% %                                   %  y1 y2 y3 ... yn];
% % end
% % % pos_array = [pos_array; ones(1,size(pos_array,2))];
% % % 2) the z depth of each points
% % z_array = zeros(1,length(pos_array));
% % for i = 1:length(z_array)
% %     temp = -sqrt(r^2 - pos_array(1,i)^2);
% %     if isreal(temp)
% %         z_array(i) = temp; %circle with r
% %     else
% %         z_array(i) = 0;
% %         disp('z is imaginary');
% %     end
% % end
% % %wheel 1
% % wheel1 = PosArray3(pos_array,z_array,1.5,2.5,pi/2);
% % % PartVisualization(wheel1_PosArray3,'wheel1');
% % Bull.PosArray3_cell{1} = wheel1;
% % 
% % wheel_x = 1.5; wheel_y = -2.5; wheel_th = -pi/4; %wheel 2
% % Bull.PosArray3_cell{2} = PosArray3(pos_array,z_array,wheel_x,wheel_y,wheel_th);
% % wheel_x = -1.5; wheel_y = -2.5; wheel_th = pi/2; %wheel 3
% % Bull.PosArray3_cell{3} = PosArray3(pos_array,z_array,wheel_x,wheel_y,wheel_th);
% % wheel_x = -1.5; wheel_y = 2.5; wheel_th = pi/2; %wheel 4
% % Bull.PosArray3_cell{4} = PosArray3(pos_array,z_array,wheel_x,wheel_y,wheel_th);
% % 
% % %for visualization : using PoseContext
% % Body = PoseContext([],Bull.get_x(),Bull.get_y(),Bull.get_z(),Bull.get_th());
% % line = [0 0;0 0;-1 1;1 1];
% % Wheel1_axis = PoseContext([],1,1,0,0);
% % Wheel1_axis.parent = Body;
% % Body.child = Wheel1_axis;
% % Wheel1 = PoseContext(line,0,0,0,0);
% % Wheel1.parent = Wheel1_axis;
% % Wheel1_axis.child = Wheel1;
% % 
% % p_world = Body.HT*Wheel1_axis.HT*Wheel1.HT*Wheel1.particles;
% % % r = 0.4;
% % % dx = 0.1;
% % % min_x = -r+dx/4;
% % % max_x = r-dx/4;
% % % x = min_x:dx/2:max_x;
% % % z = zeros(size(x));
% % % for i=1:length(x)
% % %     z(i) = -sqrt(r^2 - x(i)^2);
% % % end
% % % plot(x,z)
% 
% %%
% % dx = 0.1;
% % edges{1} = Line(0,0,1,0,dx); %x1 y1 x2 y2
% % edges{2} = Line(1,0,0,1,dx);
% % edges{3} = Line(0,1,0,0,dx);
% % 
% % % for i = 1:3
% % %     edges{i}.synchronize(dx)
% % % end
% % 
% % area_idx = get_coord_inside(edges)
% % 
% % 
% % function area_idx = get_coord_inside(lines)
% % 
% % if ~isa(lines,'cell')
% %     error('Input must be cell array of Line s');
% % end
% % dx = lines{1}.dx;
% % edges_idx = [];
% % for i = 1:3
% %     edges_idx = [edges_idx get_coord(lines{i}.points,dx)];
% % end
% % edges_idx_sorted = sortrows(edges_idx')';
% % area_idx = edges_idx; %initialize with line indices
% % for x_idx = min(edges_idx_sorted(1,:)):max(edges_idx_sorted(1,:))
% %     % from left to right
% %     temp = edges_idx_sorted(:,edges_idx_sorted(1,:)==x_idx);
% %     % points corresponding to x_idx == 1, .. 2,3.. so on
% %     % [2 2      or     [1 1 1 ... 1
% %     %  1 10]            1 2 3 ... 11]
% %     y_idx = min(temp(2,:));
% %     while y_idx < max(temp(2,:))
% %         % points between (x_idx==1,y_min) and (x_idx==1,y_max)
% %         if ~ismember(y_idx,temp(2,:))
% %             area_idx = horzcat(area_idx,[x_idx;y_idx]);
% %         end
% %         y_idx = y_idx + 1;
% %     end
% % end
% % end
% % 
% % function indices = get_coord(p,dx)
% % x_idx = round(p(1,:)/dx)+1;
% % y_idx = round(p(2,:)/dx)+1;
% % indices = [x_idx; y_idx];
% % end
% % 
% % 
% % function h = get(x)
% % h = x;
% % end
