clear all
close all
clc
wheel = matfile('wheel.mat');
area = matfile('area.mat');
normal = matfile('normal.mat');

pointList = wheel.pointList;
areaList = area.areaList;
normalList = normal.normalList;
rList = sqrt(pointList(2,:).^2+pointList(3,:).^2);

% mm/s
w = 2;

angleList = atan2(pointList(3,:),pointList(2,:))+pi/2;


vcorx = 0;
vcory = 0;
vcorz = 0;
vx = zeros([1,size(angleList,2)])+vcorx;
vz = sin(angleList).*rList.*w+vcory;
vy = cos(angleList).*rList.*w+vcorz;



figure
plot3(pointList(1,:),pointList(2,:),pointList(3,:),'ok','MarkerFaceColor',[0,0.5,0.5])
hold on
quiver3(pointList(1,:),pointList(2,:),pointList(3,:),normalList(1,:),normalList(2,:),normalList(3,:),10,'Color', [0,0.2,0.8]);
% hold on
% axisP = plot3([0,60],[0,0],[0,0]);
% axisP.LineWidth = 3;
view(45,25)
daspect([1 1 1])

% figure
% for k =1:100:size(pointList,2)
%     plot3(pointList(1,k),pointList(2,k),pointList(3,k),'ok','MarkerFaceColor',[0,0.5,0.5])
%     hold on 
%     text(pointList(1,k),pointList(2,k),pointList(3,k),string(angleList(k)*180/pi))
%     
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),vx(k),vy(k),vz(k),0.1,'r');
% end
% daspect([1 1 1])

% figure
% quiver3(pointList(1,:),pointList(2,:),pointList(3,:),vx,vy,vz,2,'r');
% view(45,25)
% daspect([1 1 1])