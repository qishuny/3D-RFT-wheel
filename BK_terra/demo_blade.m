%% -- Bulldozer Demo---
clear 
close all
n = 8; %grid per m
minx = -8;
miny = -8;
maxx = 8;
maxy = 8;
Sand = Hmap(minx,miny,maxx,maxy,n);
Sand.dT = 2;


%----------------setting up body position
bullxpos = 0;
bullypos = -6;
bullzpos = 0;
bulltheta = 0;
BodyPose = PoseContext([],bullxpos,bullypos,bullzpos,bulltheta);


%----------------setting pointcloud for blade
% set particlepoints for my wheel points wrt origin of blade pose
% 1) the x,y position of each points
x_width = 0.5;  %1/2 width of blade
xgrid = -x_width:Sand.get_dx/2:x_width; %particles' xpos
y_length = Sand.get_dx; %1/2 length of blade just enough to cover like 2~3 grids points
ygrid = -y_length:Sand.get_dx/2:y_length; %particles' ypos
xy_vec = []; 
for i = 1:length(ygrid)
    temp = [xgrid; ones(size(xgrid))*ygrid(i)];
    xy_vec = [xy_vec temp]; % [x1 x2 x3 ... xn;
                                  %  y1 y2 y3 ... yn];
end
% 2) the z position of each points
z_depth = -0.5;
z_vec = z_depth*ones(1,length(xy_vec));
particlepoints = [xy_vec;z_vec;ones(size(z_vec))];

bladexpos = 0; %always relative position of pose to parent
bladeypos = 0;
bladezpos = 0;
bladetheta = 0;
BladePose = PoseContext(particlepoints,bladexpos,bladeypos,bladezpos,bladetheta);
BladePose.parent = BodyPose;
BodyPose.child = BladePose; 

% --------------- setting lines for visualization
% blade, we need 12 lines to draw blade 
z_height = 1; % this is only for visualization
lines = cell(1,12);
w = x_width;
l = y_length;
h = z_height;
d = z_depth;
lines{1} = [-w w;  -l -l; d d; 1 1];
lines{2} = [w w;   -l -l; d h; 1 1];
lines{3} = [w -w;  -l -l; h h; 1 1];
lines{4} = [-w -w; -l -l; h d; 1 1];

lines{5} = [-w w;  l l;  d d; 1 1];
lines{6} = [w w;  l l;  d h; 1 1];
lines{7} = [w -w;  l l;  h h; 1 1];
lines{8} = [-w -w;  l l;  h d; 1 1];

lines{9} = [w w;  -l l;  d d; 1 1];
lines{10} = [w w;  -l l;  h h; 1 1];
lines{11} = [-w -w;  -l l;  d d; 1 1];
lines{12} = [-w -w;  -l l;  h h; 1 1];

BladePose.lines = lines;

%-----------initialize visualizeContext
visualizeContext = VisualizeContext(Sand.minx,Sand.miny,Sand.maxx,Sand.maxy,Sand.n);

% ----------initialize indexHandler
indexHandler = IndexHandler(Sand.get_dx,Sand.minx,Sand.miny,Sand.maxx,Sand.maxy);

% ----------initialize TrajectoryContext
trajectoryContext = TrajectoryContext();

% ----------making trajectory
initial_pos = [bullxpos bullypos bullzpos]';
initial_vel = [0 1 0]';
initial_acc = [0 0 0]';
final_pos = [bullxpos-4 bullypos+8 bullzpos]';
final_vel = [0 1 0]';
final_acc = [0 0 0]';
tfinal = 20; %10;
dT = 0.05;
% pos0 vel0 acc0, posf velf accf
[X,Y,Z,dX,dY,dZ,theta] = trajectoryContext.minJerk(initial_pos,initial_vel,initial_acc,...
                                                final_pos,final_vel,final_acc,tfinal,dT);

initial_pos = final_pos;
initial_vel = final_vel;
initial_acc = final_acc;
final_pos = initial_pos + [-2 -7 0];
final_vel = [-1 0 0]';
final_acc = [0 0 0]';
tfinal = 20; %10;
[X2,Y2,Z2,dX2,dY2,dZ2,theta2] = trajectoryContext.minJerk(initial_pos,initial_vel,initial_acc,...
                                                final_pos,final_vel,final_acc,tfinal,dT);

X = [X; X2]; Y = [Y;Y2]; Z = [Z;Z2]; dX = [dX;dX2]; dY = [dY;dY2]; dZ = [dZ;dZ2]; theta = [theta;theta2];

%-----------Resolve collision, visualize blade, the initial step
particlepointsWorld = BodyPose.HT4*BladePose.HT4*BladePose.particles;
[indices, zvec] = indexHandler.getIndex(particlepointsWorld);
velocity = [0;1];
Sand.resolve_collision(indices,zvec,velocity);
Sand.Steadystate(indices, zvec, velocity);

lines_wcoord = cell(size(lines)); %getting the world coordinate lines for visualization
for i = 1:length(BladePose.lines)
    lines_wcoord{i} = BodyPose.HT4*BladePose.HT4*BladePose.lines{i};
end
visualizeContext.visualize(Sand.matrix,lines_wcoord);

%% ----------- MAIN SIMULATION ---------------
%------------Resolve collision, visualize blade, following steps
savedataSand = cell(1,length(X));
savedatalines = cell(1,length(X));

Sand.optOn = true; 
tic
for i = 1:length(X)
    % body's trajectory
    r = sqrt(dX(i)^2+dY(i)^2);
    BodyPose.HT4 = [dY(i)/r dX(i)/r 0 X(i);
        -dX(i)/r dY(i)/r 0 Y(i);
        0 0 1 0;
        0 0 0 1];

    particlepointsWorld = BodyPose.HT4*BladePose.HT4*BladePose.particles;
    [indices, zvec] = indexHandler.getIndex(particlepointsWorld);
    velocity = [dX(i),dY(i)];
    Sand.resolve_collision(indices,zvec,velocity);
    
    Sand.Steadystate(indices, zvec, velocity); 
%     Sand.steady_state('off');
    lines_wcoord = cell(size(lines)); %getting the world coordinate lines for visualization
    for j = 1:length(BladePose.lines)
        lines_wcoord{j} = BodyPose.HT4*BladePose.HT4*BladePose.lines{j};
    end
    visualizeContext.visualize(Sand.matrix,lines_wcoord);
    savedataSand{i} = Sand.matrix;
    savedatalines{i} = lines_wcoord;
end
toc
%-----------initialize visualizeContext
% visualizeContext2 = VisualizeContext(Sand.minx,Sand.miny,Sand.maxx,Sand.maxy,Sand.n);
% 
% figure
% for i = 1:length(savedataSand)
%     visualizeContext2.visualize(savedataSand{i},savedatalines{i});
% end
