% 3D RFT Modeling on Wheels
% Qishun Yu Catherine Pavlov
% 02/21/2021

% points (3*n) data.Points
% area (1*n) data.Area
% normal vectors (3*n) data.Normals

data = matfile('data/smooth_wheel_125.mat');
% data = matfile(data/wheel_106.mat');
% data = matfile('data/plate.mat');

pointList = data.Points;
areaList = data.Area;
normalList = data.Normals;
pointSize = size(pointList, 2);


%% SET parameters
% SET slip angle
slipAngle = pi/4;
% SET velocity of the center of rotation of the body mm/s
vcenter = 10;
% SET wheel rotational speed mm/s
wr = 10;
% SET sinkage mm
sinkage = 40;
% SET radius
radius = 62.5;
% SET scaling factor
sf = 1;

%% Geometry & Velocity Calc

% angular velocity radius/s
vcorx = -vcenter * sin(slipAngle);
vcory = vcenter * cos(slipAngle);
vcorz = 0;
vcor = [vcorx; vcory; vcorz];
w = -wr / radius;

% depth = sand with respect to the center of the wheel mm
depth = -radius + sinkage;

[e1List, e2List, vList, vHoriList, v1List, v23List] = calc_velocity(normalList, pointList, w, vcor);
[phi] = calc_Phi(vHoriList, e2List);


%% Force Calc
[betaList, gammaList] = calc_BetaGamma(normalList, e2List, v23List);
[Force,netForce, idx] = calc_3D_rft(pointList, betaList, gammaList, e1List, e2List, areaList, phi, depth, sf);
Force


%% Plot stuff

% plot velocity
% figure
% quiver3(pointList(1,:),pointList(2,:),pointList(3,:),vList(1,:),vList(2,:),vList(3,:),2,'Color', [0,0.2,0.8]);
% daspect([1 1 1])


%
% figure
% for k =1:200:size(pointList,2)
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),e2List(1,k),e2List(2,k),e2List(3,k),10,'g');
%     hold on
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),v23List(1,k),v23List(2,k),v23List(3,k),0.1,'r');
%     legend('e2','v23')
%     text(pointList(1,k),pointList(2,k),pointList(3,k),string(gammaList(k)*180/pi))
%     daspect([1 1 1])
% end
%
% figure
% for k =1:200:size(pointList,2)
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),e2List(1,k),e2List(2,k),e2List(3,k),10,'g');
%     hold on
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),normalList(1,k),normalList(2,k),normalList(3,k),10,'r');
%     legend('e2','normal')
%     text(pointList(1,k),pointList(2,k),pointList(3,k),string(betaList(k)*180/pi))
%     daspect([1 1 1])
% end

%Plot selected normal vector and e1,e2

% figure
% for k =1:200:size(pointList,2)
%     plot3(pointList(1,k),pointList(2,k),pointList(3,k),'ok','MarkerFaceColor',[0,0.5,0.5])
%     hold on
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),normalList(1,k),normalList(2,k),normalList(3,k),10,'Color', [0,0.2,0.8]);
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),e1List(1,k),e1List(2,k),e1List(3,k),10,'r');
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),e2List(1,k),e2List(2,k),e2List(3,k),10,'g');
%     legend('point','normal vector','e1 axis','e2 axis')
%     daspect([1 1 1])
% end


%Plot selected velocity and v1,v23

figure
for k =1:200:size(pointList,2)
    
    quiver3(pointList(1,k),pointList(2,k),pointList(3,k),vList(1,k),vList(2,k),vList(3,k),1,'Color', [0,0.2,0.8]);
    hold on
    quiver3(pointList(1,k),pointList(2,k),pointList(3,k),v1List(1,k),v1List(2,k),v1List(3,k),1,'r');
    quiver3(pointList(1,k),pointList(2,k),pointList(3,k),v23List(1,k),v23List(2,k),v23List(3,k),1,'g');
    legend('velocity','v1','v23')
    daspect([1 1 1])
end




% plot velocity
% figure
% for k =1:150:size(forcePoints,2)
% 
%     quiver3(forcePoints(1,k),forcePoints(2,k),forcePoints(3,k),forceV1(1,k),forceV1(2,k),forceV1(3,k),0.5,'Color', [0,0.2,0.8]);
%     hold on
%     quiver3(forcePoints(1,k),forcePoints(2,k),forcePoints(3,k),forcee1(1,k),forcee1(2,k),forcee1(3,k),1,'r');
%     text(forcePoints(1,k),forcePoints(2,k),forcePoints(3,k),string(forceGamma(k)*180/pi))
%     legend('V1','e1')
% 
% end
% daspect([1 1 1])


%Plot velocity and v1 v23
% figure
% for k =1:150:size(forcePoints,2)
% 
%     quiver3(forcePoints(1,k),forcePoints(2,k),forcePoints(3,k),forceV23(1,k),forceV23(2,k),forceV23(3,k),0.5,'Color', [0,0.2,0.8]);
%     hold on
%     quiver3(forcePoints(1,k),forcePoints(2,k),forcePoints(3,k),forcee2(1,k),forcee2(2,k),forcee2(3,k),1,'r');
%     text(forcePoints(1,k),forcePoints(2,k),forcePoints(3,k),string(forceGamma(k)*180/pi))
%     legend('V23','e2')
% 
% end
% daspect([1 1 1])


% figure
% for k =1:150:size(forcePoints,2)
% 
%     quiver3(forcePoints(1,k),forcePoints(2,k),forcePoints(3,k),forceNormals(1,k),forceNormals(2,k),forceNormals(3,k),'Color', [0,0.2,0.8]);
%     hold on
%     quiver3(forcePoints(1,k),forcePoints(2,k),forcePoints(3,k),forcee2(1,k),forcee2(2,k),forcee2(3,k),1,'r');
%     text(forcePoints(1,k),forcePoints(2,k),forcePoints(3,k),string(forceBeta(k)*180/pi))
%     legend('normal','e2')
% end
% daspect([1 1 1])




% figure
% for k =1:25:size(forcePoints,2)
%     plot3(forcePoints(1,k),forcePoints(2,k),forcePoints(3,k),'ok','MarkerFaceColor',[0,0.5,0.5])
%     hold on
% %     text(forcePoints(1,k),forcePoints(2,k),forcePoints(3,k),string(forceBeta(k)*180/pi))
%     text(forcePoints(1,k),forcePoints(2,k),forcePoints(3,k),string(forceGamma(k)*180/pi))
% %     hold on
%     quiver3(forcePoints(1,k),forcePoints(2,k),forcePoints(3,k),forceV23(1,k),forceV23(2,k),forceV23(3,k),0.05,'Color', [0,0.2,0.8]);
% end
% daspect([1 1 1])

% plot alphaz
% figure
% for k = 1 : 25 : size(forcePoints, 2)
%     plot3(forcePoints(1,k),forcePoints(2,k),forcePoints(3,k),'ok','MarkerFaceColor',[0,0.5,0.5])
%     hold on 
%     text(forcePoints(1, k),forcePoints(2, k),forcePoints(3, k),string(az23(1, k)));
% end
% daspect([1 1 1]);



% plot force
figure
quiver3(pointList(1, idx),pointList(2, idx),pointList(3, idx),netForce(1,:),netForce(2,:),netForce(3,:),2,'Color', [0,0.2,0.8]);
legend('force')
hold on 
quiver3(zeros(3,1),zeros(3,1),zeros(3,1),[10;0;0],[0;10;0],[0;0;10]);
legend('unit vector')
daspect([1 1 1])
%

% figure
% for k =1:150:size(pointList,2)
%     plot3(pointList(1,k),pointList(2,k),pointList(3,k),'ok','MarkerFaceColor',[0,0.5,0.5])
%     hold on
%     text(pointList(1,k),pointList(2,k),pointList(3,k),string(phi(k)*180/pi))
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),e2List(1,k),e2List(2,k),e2List(3,k),10,'Color', [0,0.2,0.8]);
%     quiver3(pointList(1,k),pointList(2,k),pointList(3,k),vHoriList(1,k),vHoriList(2,k),vHoriList(3,k),0.1,'g');
% end
% daspect([1 1 1])

%% Functions

function [Force,netForce, idx] = calc_3D_rft(pointList, betaList, gammaList, e1List, e2List, areaList, phi, depth, sf);

idx = pointList(3,:) < depth;
forcePoints = pointList(:, idx);
numofForce = size(forcePoints, 2);

forceBeta = betaList(idx);
forceGamma = gammaList(idx);


[ay1, ~] = calc_rft_alpha(0, 0, sf);
ax23  = zeros(1,numofForce);
az23 = zeros(1,numofForce);

ay1 = zeros(1,numofForce)+ay1;

for i = 1:numofForce
    [ax23(i),az23(i)] = calc_rft_alpha(forceBeta(i),forceGamma(i), sf);
end

magF1 = -ay1 .* abs(depth - pointList(3,idx)) .* (areaList(idx) * 10^-3);
magF2 = -ax23 .* abs(depth - pointList(3,idx)) .* (areaList(idx) * 10^-3);
% F1 force in N
F1tilde = magF1 .* e1List(:,idx);
F2tilde = magF2 .* e2List(:,idx);
F2tilde(3,:) = az23 .* abs(depth - pointList(3,idx)) .* (areaList(idx) * 10^-3);


f1List = zeros(1,size(F1tilde,2));
f2List = zeros(1,size(F2tilde,2));
phiList = phi(idx);
for i = 1:size(F1tilde,2)
    [f1List(i),f2List(i)] = normalScale(phiList(i));
end

netForce = f1List .* F1tilde + f2List .* F2tilde;
Force = sum(netForce,2);
end

% find the local alphax and alphaz with give gamma and beta
% return alpha in N/(cm^3)
function [alphaX, alphaZ] = calc_rft_alpha(beta,gamma, sf)
% using discrete Fourier transform fitting function [Li et al., 2013]
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

beta = wrapToPi(beta);
gamma = wrapToPi(gamma);
if beta >= -pi && beta <= -pi/2
    beta = beta + pi;
elseif beta >= pi/2 && beta <= pi
    beta = beta - pi;
end
beta = wrapToPi(beta);
gamma = wrapToPi(gamma);
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


% Compute f1, f23 from v1 v23 v
% the sigmoid functions are from [Treers et al., 2021]
function [f1, f23] = normalScale(phi)
a1 = 0.44;
a2 = 3.62;
a3 = 1.61;
a4 = 0.41;

b1 = 1.99;
b2 = 1.61;
b3 = 0.97;
b4 = 4.31;

phi = wrapToPi(phi);

if phi>pi/2 && phi<=pi
    phi = pi-phi;
elseif phi<0 && phi>=-pi/2
    phi = -phi;
elseif phi<-pi/2 && phi >=-pi
    phi = phi+pi;
end
v1v = sin(phi);
v23v = cos(phi);

F1 = a1 * tanh(a2 * v1v - a3) + a4;
F23 = b1 * atanh(b2 * v23v - b3) + b4;
F1max =a1 * tanh(a2 * sin(pi/2) - a3) + a4;
F23max = b1 * atanh(b2 * cos(0) - b3) + b4;

f1 = F1 / F1max;
f23 = F23 / F23max;
end


function angles = calc_Angles(v1, v2)
dotprd =v1(1,:).*v2(1,:)+v1(2,:).*v2(2,:)+v1(3,:).*v2(3,:);
timeprd = sqrt(v1(1,:).^2+v1(2,:).^2+v1(3,:).^2).*sqrt(v2(1,:).^2+v2(2,:).^2+v2(3,:).^2);

% not stable, may output imaginary value when angle is around pi. 
angles = acos((dotprd)./(timeprd));
end


function [e1List, e2List, vList, vHoriList, v1List, v23List] = calc_velocity(normalList, pointList, w, vcor)
numofNormal = size(normalList, 2);
%e2 unit axis
e2List = [normalList(1, :);
    normalList(2, :);
    zeros(1, numofNormal)];
idxe2 = (e2List(1, :) == 0 & e2List(2, :) == 0);
e2List(2,:) = 1 .* idxe2 + e2List(2,:) .* (~idxe2);
magne2 = sqrt(e2List(1, :) .^2 + e2List(2, :) .^2 + e2List(3, :) .^2);
e2List = [e2List(1, :) ./ magne2;
    e2List(2, :) ./ magne2;
    e2List(3, :) ./ magne2;];

% element radius to the center information (1*n)
rList = sqrt(pointList(2,:).^2+pointList(3,:).^2);   
angleList = atan2(pointList(3, :), pointList(2, :)) + pi/2;

vx = zeros([1, size(angleList, 2)]) + vcor(1);
vy = cos(angleList) .* rList .* w + vcor(2);
vz = sin(angleList) .* rList .* w + vcor(3);

% velocity
vList = [vx; 
    vy; 
    vz];

vHoriList = [vx; vy; zeros(1, size(angleList, 2))];
idxV = (vHoriList(1, :) == 0 & vHoriList(2, :) == 0);
vHoriList(2, :) = 1 .* idxV + vHoriList(2, :) .* (~idxV);
vHorie2List = dot(vList, e2List) ./ 1 .* e2List;

% e1 unit axis
e1List = vHoriList - vHorie2List;
magne1 = sqrt(e1List(1, :) .^2 + e1List(2, :) .^2 + e1List(3, :) .^2);
e1List = [e1List(1, :) ./ magne1;
    e1List(2, :) ./ magne1;
    e1List(3, :) ./ magne1;];

% v1 and v23
v1List = dot(vList, e1List) ./ 1 .* e1List;
v23List = vList - v1List;

end


function [phi] = calc_Phi(vHoriList, e2List)
phi = calc_Angles(vHoriList, e2List);
phi = wrapToPi(phi);
end


function [betaList, gammaList] = calc_BetaGamma(normalList, e2List, v23List)
betaList = calc_Angles(normalList, e2List);
idxBeta = normalList(3,:) < 0;
betaList(idxBeta) = - betaList(idxBeta);
betaList = betaList + pi/2;
betaList = wrapToPi(betaList);

gammaList = calc_Angles(v23List, e2List);
idxGamma = v23List(3,:) > 0;
gammaList(idxGamma) = - gammaList(idxGamma);
gammaList = wrapToPi(gammaList);

end


%UNUSED FUNCTIONS

% [e1List, e2List] = calc_e1e2(normalList);
% [vList, vHoriList, v1List, v23List] = calc_Vel(pointList, w, radius, e1List, vcor);

function [e1List, e2List] = calc_e1e2(normalList)
numofNormal = size(normalList, 2);
%e2 unit axis
e2List = [normalList(1, :);
    normalList(2, :);
    zeros(1, numofNormal)];
idxe2 = (e2List(1, :) == 0 & e2List(2, :) == 0);
e2List(2,:) = 1 .* idxe2 + e2List(2,:) .* (~idxe2);

magne2 = sqrt(e2List(1, :) .^2 + e2List(2, :) .^2 + e2List(3, :) .^2);
e2List = [e2List(1, :) ./ magne2;
    e2List(2, :) ./ magne2;
    e2List(3, :) ./ magne2;];
%e1 unit axis
e1List = [e2List(1,:).*cos(-pi/2)+e2List(2,:).*-sin(-pi/2);
    e2List(1,:).*sin(-pi/2)+e2List(2,:).*cos(-pi/2);
    e2List(3,:)];

end
function [vList, vHoriList, v1List, v23List] = calc_Vel(pointList, w, radius,e1List,vcor)
% element radius to the center information (1*n)
rList = sqrt(pointList(2,:).^2+pointList(3,:).^2);

% element tangent angle (1*n)
angleList = atan2(pointList(3,:),pointList(2,:))+pi/2;

vx = zeros([1,size(angleList,2)])+vcor(1);
vy = cos(angleList).*rList.*w+vcor(2);
vz = sin(angleList).*rList.*w+vcor(3);

vList = [vx;vy;vz];

vHoriList = [vx;vy;zeros(1,size(angleList,2))];

idxV = (vHoriList(1,:)==0 & vHoriList(2,:)==0);
vHoriList(2,:) = 1.*idxV+vHoriList(2,:).*(~idxV);
v1List = dot(vList,e1List)./1.*e1List;
v23List = vList-v1List;
end
