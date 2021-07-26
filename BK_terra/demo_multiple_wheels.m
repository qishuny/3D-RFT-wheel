%% -- Multiple Wheel Demo---
clear 
close all

% parameter for wheel 
% diameter: 0.096m
% width of wheel: 0.048m
% angle of repose 29 degrees

n = 200; %grid per m; make it at most 0.02
minx = -0.4;
miny = -0.4;
maxx = 0.4;
maxy = 0.4;
Sand = Hmap(minx,miny,maxx,maxy,n);


%----------------setting up body position
bodyxpos = 0;
bodyypos = miny + abs(miny/5); %miny + abs(miny/5);
bodyzpos = 0.096/2-0.012;
bodytheta = 0;
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

% ---- first wheel
% join two poses: body and wheel
wheelxpos = -0.2; %always relative position of pose to parent
wheelypos = 0;
wheelzpos = 0;
wheeltheta = 45*pi/180;
WheelPose = PoseContext(particlepoints,wheelxpos,wheelypos,wheelzpos,wheeltheta);
WheelPose.parent = BodyPose;
BodyPose.child = WheelPose; 

% ---- second wheel
wheelxpos = 0.2; %always relative position of pose to parent
wheelypos = 0;
wheelzpos = 0;
wheeltheta = 90*pi/180;
WheelPose2 = PoseContext(particlepoints,wheelxpos,wheelypos,wheelzpos,wheeltheta);
WheelPose2.parent = BodyPose;
BodyPose.child = WheelPose2; 


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
WheelPose2.lines = lines;


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
final_pos = [bodyxpos bodyypos+(maxy-initial_pos(2))-maxy/3 bodyzpos]'; %bodyypos+(maxy-initial_pos(2))-maxy/3
final_vel = [0 initial_vel(2) 0]';
final_acc = [0 0 0]';
tfinal = norm(final_pos - initial_pos)/norm(initial_vel); %10;
dT = 1;
% gives me vector of X,Y,Z for trajectory
[X,Y,Z,dX,dY,dZ,theta] = trajectoryContext.minJerk(initial_pos,initial_vel,initial_acc,final_pos,final_vel,final_acc,tfinal,dT);

%-----------Resolve collision, visualize blade, the initial step
particlepointsWorld = BodyPose.HT4*WheelPose.HT4*WheelPose.particles;
temp = BodyPose.HT4*WheelPose2.HT4*WheelPose2.particles;
particlepointsWorld = [particlepointsWorld temp];
[indices, zvec] = indexHandler.getIndex(particlepointsWorld);
velocity = [0;1];
Sand.resolve_collision(indices,zvec,velocity);
Sand.Steadystate(indices, zvec, velocity); % Steadystate wayyyy too slow

% for visualization 
lines_wcoord = cell(size(lines)); %getting the world coordinate lines for visualization
for i = 1:length(WheelPose.lines)
    lines_wcoord{i} = BodyPose.HT4*WheelPose.HT4*WheelPose.lines{i};
end
for i = 1:length(WheelPose2.lines)
    lines_wcoord{end+1} = BodyPose.HT4*WheelPose2.HT4*WheelPose2.lines{i};
end

visualizeContext.visualize(Sand.matrix,lines_wcoord);
Sand.Steadystate(indices, zvec, velocity); 
visualizeContext.visualize(Sand.matrix,lines_wcoord);


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
    temp = BodyPose.HT4*WheelPose2.HT4*WheelPose2.particles;
    particlepointsWorld = [particlepointsWorld temp];
    
    [indices, zvec] = indexHandler.getIndex(particlepointsWorld);
    velocity = [dX(i),dY(i)];

    Sand.Steadystate(indices,zvec,velocity);
    Sand.resolve_collision(indices,zvec,velocity);


    lines_wcoord = cell(size(lines)); %getting the world coordinate lines for visualization
    for j = 1:length(WheelPose.lines)
        lines_wcoord{j} = BodyPose.HT4*WheelPose.HT4*WheelPose.lines{j};
    end
    for i = 1:length(WheelPose2.lines)
        lines_wcoord{end+1} = BodyPose.HT4*WheelPose2.HT4*WheelPose2.lines{i};
    end

    visualizeContext.visualize(Sand.matrix,lines_wcoord);
    savedataSand{i} = Sand.matrix;
    savedatalines{i} = lines_wcoord;
end
