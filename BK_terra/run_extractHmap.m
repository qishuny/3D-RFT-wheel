wheeldata = matfile('smooth_wheel_125.mat');
pointList = wheeldata.Points;
normalList = wheeldata.Normals;



pointList(1,:) = pointList(1,:) - 30;
pointListOriginal = pointList;
% slipAngle in degree
slipAngle = 90;

wheelDiameter = 0.125; %m
wheelWidth = 0.06; %m
depth = 0.04; %m

%grid per m
n = 400;

gridsize = 1/n;

[pointList, normalList] = rotateZ(pointList, normalList, -slipAngle * pi / 180);


tic
[sandHmap, wheelPos] = extractHmap((90 - slipAngle), wheelDiameter, wheelWidth, depth, n);
sandHmap = sandHmap';
toc

%% change unit to mm
sandHmap = sandHmap * 1000;
wheelPos = wheelPos * 1000;
wheelDiameter = wheelDiameter * 1000;
wheelWidth = wheelWidth * 1000;
depth = depth * 1000; 
gap = gridsize * 1000;
%% plot simluation result: sand height map & wheel position
close all

[X, Y] = meshgrid(-400:gap:400, -400:gap:400);


figure
s = surf(X, Y, sandHmap);
s.EdgeColor = 'none';
hold on
scatter3(wheelPos(1), wheelPos(2), wheelPos(3),'r')
axis on
xlabel('x')
ylabel('y')
axis equal




%% line up the wheel and the sand height map

lim = sqrt((wheelDiameter / 2) ^ 2 + (wheelWidth / 2) ^ 2) + 10;
idxX = X(1, :) >(wheelPos(1) - lim) & X(1, :) < (wheelPos(1) + lim);
idxY = Y(:, 1) >(wheelPos(2) - lim) & Y(:, 1) < (wheelPos(2) + lim);

Y = Y - wheelPos(2);
sandHmap = sandHmap - wheelPos(3);

wheelPos(3) = 0;

Xtrimed = X(idxY, idxX);
Ytrimed = Y(idxY, idxX);
SandHmapOriginal = sandHmap(idxY, idxX);


figure

plot3(pointList(1,:), pointList(2,:), pointList(3,:), '.', ...
    'Color', [0.0,0.0,0.0], 'MarkerSize', 1)
hold on
surf(Xtrimed, Ytrimed, SandHmapOriginal, 'FaceAlpha', 0.5)

axis on
xlabel('x')
ylabel('y')
axis equal

%% 


XY = [Xtrimed(:) Ytrimed(:)];  
theta = slipAngle;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
rotXY = XY * R';
Xqr = reshape(rotXY(:,1), size(Xtrimed,1), []);
Yqr = reshape(rotXY(:,2), size(Ytrimed,1), []);






idxWheelX = Xqr(:, :) <=  wheelWidth / 2 + 2 & Xqr(:, :) >= - wheelWidth/ 2 - 2;
idxWheelY = Yqr(:, :) <= wheelDiameter/ 2  + 2& Yqr(:, :) >= - wheelDiameter/ 2 -2;

idx = idxWheelX & idxWheelY;
SandHmapnew = SandHmapOriginal;

SandHmapnew(idx) = max(SandHmapnew(idx), -depth);
% SandHmapnew(idx) = -depth;

figure
surf(Xqr, Yqr, SandHmapnew)

[Points_inflated] = inflate(pointListOriginal, wheelDiameter, wheelWidth, gap);
[Points_inflated, normalList] = rotateZ(Points_inflated, normalList, -slipAngle * pi/180);

spz = interp2(Xtrimed, Ytrimed, SandHmapnew, Points_inflated(1,:)', Points_inflated(2,:)');

mask = Points_inflated(3,:) > spz';
figure
plot3(Points_inflated(1,mask), Points_inflated(2,mask), Points_inflated(3,mask),'.','Color',[0.8,0.8,0.8],'MarkerSize',1);
hold on
plot3(Points_inflated(1,~mask), Points_inflated(2,~mask), Points_inflated(3,~mask),'.','Color','r','MarkerSize',1);
s = surf(Xtrimed, Ytrimed, SandHmapnew, 'FaceAlpha', 0.5);
s.EdgeColor = 'none';
axis equal


figure
plot3(Points_inflated(1,mask), Points_inflated(2,mask), Points_inflated(3,mask),'.','Color',[0.8,0.8,0.8],'MarkerSize',1);
hold on
plot3(Points_inflated(1,~mask), Points_inflated(2,~mask), Points_inflated(3,~mask),'.','Color','r','MarkerSize',1);
axis equal

function [Points, Normals] = rotateZ(Points, Normals,thetaz)

Rz = [cos(thetaz), -sin(thetaz), 0;
    sin(thetaz), cos(thetaz), 0;
    0, 0, 1;];
Points = Rz * Points;
Normals = Rz * Normals;
end

function [Points_inflated] = inflate(Points, wheelDiameter, wheelWidth, gap)
scaleY = (wheelDiameter / 2 + gap) / (wheelDiameter /2);
scaleX = (wheelWidth / 2 + gap) / (wheelWidth /2);

Points_inflated(1,:) = Points(1,:) .* scaleX;
Points_inflated(2,:) = Points(2,:) .* scaleY;
Points_inflated(3,:) = Points(3,:);
end 