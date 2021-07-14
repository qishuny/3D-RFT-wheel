close all

wheeldata = matfile('data/grousered_wheel_125.mat');
pointList = wheeldata.Points;
pointList(1,:) = pointList(1,:) - 30;
areaList = wheeldata.Area;
normalList = wheeldata.Normals;
sinkage = 20;
slipAngle = pi/4;
radius = 62.5;
depth = -radius + sinkage;
[pointList, normalList] = rotateZ(pointList, normalList, -slipAngle);

figure 
idx = pointList(3,:) > depth;
pointList1 = pointList(:,idx);
pointList2 = pointList(:,~idx);
plot3(pointList1(1,:),pointList1(2,:),pointList1(3,:),'.','Color',[0.9,0.9,0.9],'MarkerSize',1)

hold on
plot3(pointList2(1,:),pointList2(2,:),pointList2(3,:),'.','Color',[0.0,0.0,0.0],'MarkerSize',1)


[X,Y] = ndgrid(-80:2:80,-80:2:80);
Z = 0 .* X .* Y;
Z(:,:) = depth;
plot3(X,Y,Z,'Color', 'y');
% x = [-70 -70 70 70 -70];
% y = [-70 70 70 -70 -70];
% z = [depth depth depth depth depth];
% plot3(x, y, z,'y')
% patch(x,y,z,'yellow')

h = arrow([0, 0, 0],[0, 80, 0]);
view(-100,15)

% rectangle('Position',[-100,-100,200,200],'FaceColor','y','edgecolor','y')
xlim([-80 80])
ylim([-80 80])
xlabel('X')
ylabel('Y')
zlabel('Z')
daspect([1 1 1])
function [Points, Normals] = rotateZ(Points, Normals,thetaz)

Rz = [cos(thetaz), -sin(thetaz), 0;
    sin(thetaz),cos(thetaz), 0;
    0,0, 1;];
Points = Rz * Points;
Normals = Rz * Normals;
end