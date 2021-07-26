%% -- One Wheel Demo---
function [sandHmap, wheelPos] = extractHmap(slipAngle, wheelDiameter, wheelWidth, depth, n);


% parameter for wheel 
% diameter: 0.096m
% width of wheel: 0.048m
% angle of repose 29 degrees

%grid per m


%four corners of height map
minx = -0.4; miny = -0.4;
maxx = 0.4; maxy = 0.4;

% Hmap contains height map matrix,
% updates height map until angle of repose 
Sand = Hmap(minx,miny,maxx,maxy,n);
Sand.optOn = true; 

% VisualizeContext visualizes height map and vehicles (PoseContext)
visualizeContext = VisualizeContext(Sand.minx,Sand.miny,Sand.maxx,Sand.maxy,Sand.n);

% IndexHandler x,y,z coordinates from PoseContext, Hmap, and user input
%and translates it to Sand's matrix index. This is also separated because
%index translation is used both by Vehicle and Hmap
indexHandler = IndexHandler(Sand.get_dx,Sand.minx,Sand.miny,Sand.maxx,Sand.maxy);

% TrajectoryContext contains some pre-made trajectories that I can use
% easily, I just need to tell him start and end point
trajectoryContext = TrajectoryContext();

%----------------setting up the wheel (PoseContext)
bodyxpos = 0;
bodyypos = miny + abs(miny/5); 
% 0.012m depth
bodyzpos = wheelDiameter/2 - depth;
% relative orientation of the wheel 
bodytheta = pi/4;

% PoseContext allows me to add multiple objects to the Vehicle
%and rotate Wheel without rotating the whole body. PoseContext organize
%the hiearchy of parent and child poses 
% PoseContext(particles, xpos, ypos, zpos, theta) x,y,z,th relative to its
% parent. BodyPose will have one or more children poses
BodyPose = PoseContext([],bodyxpos,bodyypos,bodyzpos,bodytheta);


%----------------setting pointcloud for wheel (collision)
% set particlepoints for my wheel wrt origin of Vehicle
% 1) the x,y position of each points

%half the radius of wheel
x_width = wheelDiameter/2;  
% wheel x points
xgrid = -x_width:Sand.get_dx/2:x_width; 

%wheel's fatness
y_length = wheelWidth/2; 
% wheel y points
ygrid = -y_length:Sand.get_dx/2:y_length; 

xy_vec = []; 
for i = 1:length(ygrid)
    temp = [xgrid; ones(size(xgrid))*ygrid(i)];
    xy_vec = [xy_vec temp]; %#ok<AGROW> 
    %xy_vec = [x1 x2 x3 ... xn;
    %          y1 y2 y3 ... yn];
end

% wheel z points
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

% join two poses: Vehicle and Wheel
wheelxpos = 0; 
wheelypos = 0;
wheelzpos = 0;
wheeltheta = slipAngle*pi/180;
%Wheel's PoseContext
WheelPose = PoseContext(particlepoints,wheelxpos,wheelypos,wheelzpos,wheeltheta);
%attach two poses
WheelPose.parent = BodyPose;
BodyPose.child = WheelPose; 


% --------------- setting lines for visualization (visualization)
% tire line
t = 0:pi/10:2*pi;
st = r*sin(t);
ct = r*cos(t);
uno = ones(size(t));

%circle line 1
x_points = ct;
y_points = y_length*uno;
z_points = st;
lines{1} = [x_points;y_points;z_points;uno];
%circle line 2
y_points = -y_length*uno;
lines{2} = [x_points;y_points;z_points;uno];
%tire lines
for i = 1:length(t)
    lines{2+i} = [ct(i), ct(i);-y_length,y_length;st(i),st(i);1 1]; 
end
%WheelPose now contains particlespoints for collision and 
%lines for visualization
WheelPose.lines = lines;

%line info for visualization which will be used by VisualizationContext
lines_wcoord = cell(size(lines));
for i = 1:length(WheelPose.lines)
    lines_wcoord{i} = BodyPose.HT4*WheelPose.HT4*WheelPose.lines{i};
end

% ----------making trajectory
initial_pos = [bodyxpos bodyypos bodyzpos]';
% velocity is half of resolution so that there is enough trajectory points 
initial_vel = [0 Sand.get_dx/2 0]';
initial_acc = [0 0 0]';
final_pos = [bodyxpos bodyypos+(maxy-initial_pos(2))-maxy bodyzpos]'; %bodyypos+(maxy-initial_pos(2))-maxy/3
final_vel = [0 initial_vel(2) 0]';
final_acc = [0 0 0]';
tfinal = norm(final_pos - initial_pos)/norm(initial_vel); %10;
% dT does not really have a meaning here. 
dT = 1;
% gives me vector of X,Y,Z for trajectory
[X,Y,Z,dX,dY,dZ,theta] = trajectoryContext.minJerk(initial_pos,initial_vel,initial_acc,final_pos,final_vel,final_acc,tfinal,dT);


%% Simulation Begins Here

%-----------Resolve collision, visualize blade, the initial step
particlepointsWorld = BodyPose.HT4*WheelPose.HT4*WheelPose.particles;
[indices, zvec] = indexHandler.getIndex(particlepointsWorld);
velocity = [0;1];

sprintf('Volume of sand before interaction is %d',sum(sum(Sand.matrix)))
%remove sand from where wheel is to along the velocity vector
Sand.resolve_collision(indices,zvec,velocity);

%after removing sand, update Sand enough times to ensure steady state. Turn off
%visualization off while doing this: this is artifacts of my old code, now
%visualization is handled by VisualizationContext only 
% Sand.steady_state('off'); 

%or update one step multiple times (this is essentially what is happening
%in Sand.steady_state
Sand.update_onestep;
for i=1:140
    Sand.update_onestep;
    %for visualizing dynamics of Sand
%     visualizeContext.visualize(Sand.matrix,lines_wcoord);
end

% visualizeContext.visualize(Sand.matrix,lines_wcoord);

%------------Resolve collision, visualize blade, following steps
%some cell for saving the process
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
%     visualizeContext.visualize(Sand.matrix,lines_wcoord);
    savedataSand{i} = Sand.matrix;
    savedatalines{i} = lines_wcoord;
end
sprintf('Volume of sand after interaction is %d',sum(sum(Sand.matrix)))

sandHmap = Sand.matrix;
wheelPos = final_pos;

% sampling part of the track for Cat

% points_from = indexHandler.getIndex([-0.2;0;0]);
% points_to = indexHandler.getIndex([0.2;0;0]);
% 
% sliceofpoints = Sand.matrix(points_from(1):points_to(1),points_from(2):points_to(2));
% xpoints = -0.2:Sand.get_dx:0.2;
% 
% 
% figure
% plot(xpoints,sliceofpoints)
% axis equal



end
