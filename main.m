clear all
close all
clc

% read element information generated from meshstuff.m
% element values are stored in .mat files
wheel = matfile('wheel.mat');
area = matfile('area.mat');
normal = matfile('normal.mat');

% WILL CHANGE WITH ROTATION
% element points information (3*n) [x;y;z]
pointList = wheel.pointList;

% WILL NOT CHANGE WITH ROTATION
% element areas information (1*n)
areaList = area.areaList;
% element radius to the center information (1*n)
rList = sqrt(pointList(2,:).^2+pointList(3,:).^2);
numofElements = size(pointList,2);

% WILL CHANGE WITH ROTATION
% element tangent angle (1*n)
angleList = atan2(pointList(3,:),pointList(2,:))+pi/2;
% element normal vectors information (3*n)
normalList = normal.normalList;
numofNormal = size(normalList,2);
%e2 unit axis
e2List = [normalList(1,:);
    normalList(2,:);
    zeros(1,numofNormal)];
magn = sqrt(e2List(1,:).^2+e2List(2,:).^2+e2List(3,:).^2);
e2List = [e2List(1,:)./magn;
    e2List(2,:)./magn;
    e2List(3,:)./magn;];
%e1 unit axis
e1List = [e2List(1,:).*cos(-pi/2)+e2List(2,:).*-sin(-pi/2);
    e2List(1,:).*sin(-pi/2)+e2List(2,:).*cos(-pi/2);
    e2List(3,:)];


% angular speed mm/s
w = 2;
% the velocity of the center of the wheel mm/s
vcorx = 200;
vcory = 0;
vcorz = 0;
vx = zeros([1,size(angleList,2)])+vcorx;
vz = sin(angleList).*rList.*w+vcory;
vy = cos(angleList).*rList.*w+vcorz;

vList = [vx;vy;vz];

v1List = dot(vList,e1List)./1.*e1List;
v23List = vList-v1List;
% add pi/2
nAngleList = atan2(normalList(3,:),sqrt(normalList(2,:).^2+normalList(1,:).^2));
betaList = nAngleList+pi/2;
betaList = wrapToPi(betaList);
gammaList = atan2(v23List(3,:),sqrt(v23List(2,:).^2+v23List(1,:).^2));
gammaList = wrapToPi(gammaList);



[ay1,az1] = alpha_func(0,0);
ax23  = zeros(1,numofNormal);
az23 = zeros(1,numofNormal);


F1 = zeros(1,numofNormal)+ay1;
F23 = zeros(1,numofNormal);

for i = 1:numofNormal
    [ax23(i),az23(i)] = alpha_func(betaList(i),gammaList(i));
end


% figure
% for k =1:150:size(pointList,2)
%     plot3(pointList(1,k),pointList(2,k),pointList(3,k),'ok','MarkerFaceColor',[0,0.5,0.5])
%     hold on 
% %     text(pointList(1,k),pointList(2,k),pointList(3,k),string(betalList(k)*180/pi))
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),normalList(1,k),normalList(2,k),normalList(3,k),10,'Color', [0,0.2,0.8]);
%     
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),e1List(1,k),e1List(2,k),e1List(3,k),10,'r');
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),e2List(1,k),e2List(2,k),e2List(3,k),10,'g');
%     legend('point','normal vector','e2 axis','e1 axis')
% end
% daspect([1 1 1])
% 
% figure
% for k =1:150:size(pointList,2)
%     plot3(pointList(1,k),pointList(2,k),pointList(3,k),'ok','MarkerFaceColor',[0,0.5,0.5])
%     hold on 
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),vList(1,k),vList(2,k),vList(3,k),0.05,'Color', [0,0.2,0.8]);
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),v1List(1,k),v1List(2,k),v1List(3,k),0.05,'r');
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),v23List(1,k),v23List(2,k),v23List(3,k),0.05,'g');
%     legend('point','velocity','v1','v23')
% end
% daspect([1 1 1])
% 
% figure
% for k =1:150:size(pointList,2)
%     plot3(pointList(1,k),pointList(2,k),pointList(3,k),'ok','MarkerFaceColor',[0,0.5,0.5])
%     hold on 
%     text(pointList(1,k),pointList(2,k),pointList(3,k),string(nAngleList(k)*180/pi))
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),normalList(1,k),normalList(2,k),normalList(3,k),10,'Color', [0,0.2,0.8]);
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),e2List(1,k),e2List(2,k),e2List(3,k),10,'g');
% end
% daspect([1 1 1])

function [alphaX, alphaZ] = alpha_func(beta,gamma)
    % using discrete Fourier transform fitting function
    % beta [-pi,pi]
    % gamma [-pi,pi]
    % Fourier coefficients M
    % granular medium: generic coefficient
    % define (-1 as M1 for simplicity)
    A00 = 0.206;
    A10 = 0.169;
    B11 = 0.212;
    B01 = 0.358;
    BM11 = 0.055;
    C11 = -0.124;
    C01 = 0.253;
    CM11 = 0.007;
    D10 = 0.088;
    M = [A00, A10, B11, B01, BM11, C11, C01, CM11, D10];
    % scaling factor
    sf = 1;

    if beta >= -pi && beta <= -pi/2
        beta = beta + pi;
    elseif beta >= pi/2 && beta <= pi
        beta = beta - pi;
    end

    if gamma >= -pi && gamma <= -pi/2
        alphaZ = sf*(M(1)*cos(0)+M(2)*cos(2*(-beta))+M(3)*sin(2*(-beta)+(-pi-gamma))...
            +M(4)*sin((-pi-gamma))+M(5)*sin((-2*(-beta))+(-pi-gamma)));
        alphaX = -sf*(M(6)*cos(2*(-beta)+(-pi-gamma))+M(7)*cos((-pi-gamma))...
            +M(8)*sin(-2*(-beta)+(-pi-gamma))+M(9)*sin(2*(-beta)));
    elseif gamma >= pi/2 && gamma <= pi
        alphaZ = sf*(M(1)*cos(0)+M(2)*cos(2*(-beta))+M(3)*sin(2*(-beta)+(pi-gamma))...
            +M(4)*sin((pi-gamma))+M(5)*sin((-2*(-beta))+(pi-gamma)));
        alphaX = -sf*(M(6)*cos(2*(-beta)+(pi-gamma))+M(7)*cos((pi-gamma))...
            +M(8)*sin(-2*(-beta)+(pi-gamma))+M(9)*sin(2*(-beta)));
    else
        alphaZ = sf*(M(1)*cos(0)+M(2)*cos(2*beta)+M(3)*sin(2*beta+gamma)...
            +M(4)*sin(gamma)+M(5)*sin((-2*beta)+gamma));
        alphaX = sf*(M(6)*cos(2*beta+gamma)+M(7)*cos(gamma)...
            +M(8)*sin(-2*beta+gamma)+M(9)*sin(2*beta));
    end

end