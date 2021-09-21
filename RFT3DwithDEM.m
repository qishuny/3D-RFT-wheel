% 3D RFT Modeling on Wheels
% Qishun Yu Catherine Pavlov
% 02/21/2021

% global variables for calc_alpha and calc normal_scale
global M a1 a2 a3 a4 b1 b2 b3 b4
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

a1 = 0.44;
a2 = 3.62;
a3 = 1.61;
a4 = 0.41;
b1 = 1.99;
b2 = 1.61;
b3 = 0.97;
b4 = 4.31;

% load experiment data
load('output/all_smooth_data_2.mat')

% load wheel point data

wheeldata = matfile('data/smooth_wheel_125.mat');
% wheeldata = matfile('data/grousered_wheel_125.mat');
% wheeldata = matfile('data/plate.mat');

pointList = wheeldata.Points;
areaList = wheeldata.Area;
normalList = wheeldata.Normals;

% 1 for plot 0 for not plot
plotForce = 1;
plotVelocity = 0;
plotGeometry = 0;

% 1 for run all data 
runData_toggle = 1;

%% SET parameters
% SET slip angle
slipAngle = pi/4;
% SET velocity of the center of rotation of the body mm/s
vcenter = 10;
% SET wheel rotational speed mm/s
wr = 0.001;
% SET sinkage mm
sinkage = 54;
% SET radius
radius = 62.5;

% SET scaling factor
sf = 1;

%% Single RFT calc

% Geometry & Velocity Calc
tic
[e1List, e2List, vList, vHoriList, v1List, v23List, phi] = calc_velocity(normalList, pointList, wr, vcenter, radius, slipAngle);
% Force Calc
[betaList, gammaList] = calc_BetaGamma(normalList, e2List, v23List);
[Force, netForce, idx] = calc_3D_rft(pointList, betaList, gammaList, e1List, e2List, areaList, phi, sinkage, radius, sf, slipAngle);
toc

% transfer to experiment result frame
Ftractive = Force(2);
Fsidewall = -Force(1);
Fload = Force(3);

%% run all slip conditions

if runData_toggle == 1
    
    h = waitbar(0,'Initializing waitbar...');
    tic
    runData(all_results, pointList, normalList, areaList, vcenter, radius, sf, h);
    toc
end


%% plot force
if plotForce == 1
    figure
    quiver3(pointList(1, idx),pointList(2, idx),pointList(3, idx),netForce(1,:),netForce(2,:),netForce(3,:),2,'Color', [0,0.2,0.8]);

    hold on 

    legend('force');
    daspect([1 1 1])
    view(-55,15)
    axis off
end

%% plot velocity
gapSize = 20;
if plotVelocity == 1  
    %Plot selected velocity and v1,v23
    figure
    for k =1:gapSize:size(pointList,2)

        quiver3(pointList(1,k),pointList(2,k),pointList(3,k),vList(1,k),vList(2,k),vList(3,k),1,'Color', [0,0.2,0.8]);
        hold on
        quiver3(pointList(1,k),pointList(2,k),pointList(3,k),v1List(1,k),v1List(2,k),v1List(3,k),1,'r');
        quiver3(pointList(1,k),pointList(2,k),pointList(3,k),v23List(1,k),v23List(2,k),v23List(3,k),1,'g');
        legend('velocity','v1','v23')
        daspect([1 1 1])
    end

    % Plot v23 and e2
    figure
    for k =1:gapSize:size(pointList,2)
        quiver3(pointList(1,k),pointList(2,k),pointList(3,k),e2List(1,k),e2List(2,k),e2List(3,k),'g');
        hold on
        quiver3(pointList(1,k),pointList(2,k),pointList(3,k),v23List(1,k),v23List(2,k),v23List(3,k),'r');
        legend('e2','v23')
    %     text(pointList(1,k),pointList(2,k),pointList(3,k),string(gammaList(k)*180/pi))
        daspect([1 1 1])
    end

    %Plot v1 and e1
    figure
    for k =1:gapSize:size(pointList,2)
        quiver3(pointList(1,k),pointList(2,k),pointList(3,k),e1List(1,k),e1List(2,k),e1List(3,k),'g');
        hold on
        quiver3(pointList(1,k),pointList(2,k),pointList(3,k),v1List(1,k),v1List(2,k),v1List(3,k),'r');
        legend('e1','v1')
        daspect([1 1 1])
    end
end

%% Plot geometry
if plotGeometry == 1  
    %Plot selected normal vector and e1,e2
    figure
    plot3(pointList(1,:),pointList(2,:),pointList(3,:),'ok','MarkerFaceColor',[0,0.5,0.5])
    hold on
    quiver3(pointList(1,:),pointList(2,:),pointList(3,:),normalList(1,:),normalList(2,:),normalList(3,:),10,'Color', [0,0.2,0.8]);
    quiver3(pointList(1,:),pointList(2,:),pointList(3,:),e1List(1,:),e1List(2,:),e1List(3,:),10,'r');
    quiver3(pointList(1,:),pointList(2,:),pointList(3,:),e2List(1,:),e2List(2,:),e2List(3,:),10,'g');
    legend('point','normal vector','e1 axis','e2 axis')
    daspect([1 1 1])
end

function runData(all_results, pointList, normalList, areaList, vcenter, radius, sf, h)
for i=1:length(all_results)
    result = all_results(i);
    wr = result.Vry;
    if wr == 0
        wr = 0.00001;
    end  
    sinkage = abs(result.avg_Z);
%     sinkage = 20;
    slipAngle = result.beta * pi / 180;
    

    % Geometry & Velocity Calc
    [e1List, e2List, ~, ~, ~, v23List, phi] = ...
        calc_velocity(normalList, pointList, wr, vcenter, radius, slipAngle);

    % Force Calc
    [betaList, gammaList] = calc_BetaGamma(normalList, e2List, v23List);
    [Force, ~, ~] = ...
        calc_3D_rft(pointList, betaList, gammaList, e1List, e2List, areaList, phi, sinkage, radius, sf, slipAngle);

    Fsidewall = -Force(1);
    Ftractive = Force(2);
    Fload = Force(3);
    RFToutput(i) = struct('ForceX', Ftractive, ...
        'ForceY', Fsidewall , ...
        'ForceZ', Fload, ...
        'wr', result.Vry, ...
        'depth', result.avg_Z, ...
        'beta', result.beta, ...
        'slip', result.slip); 
    
    waitbar(i/length(all_results), h, 'In progress...',sprintf('%12.9f',i))
end
waitbar(1,h,'Completed.');
disp("Done.");

close(h);
save('output/RFTDEMoutput.mat','RFToutput');


end

function [Force, netForce, idx] = calc_3D_rft(pointList, betaList, gammaList, e1List, e2List, areaList, phi, sinkage, radius, sf, slipAngle);
% depth = sand with respect to the center of the wheel mm
depth = -radius + sinkage;

% find points below the surface of the soil

[idx, depthList] = run_extractHmap(pointList, slipAngle * 180 / pi, abs(sinkage / 1000));

sum(idx)
% idx = pointList(3,:) < depth;
% depthList = abs(depth - pointList(3,idx));
forcePoints = pointList(:, idx);
numofForce = size(forcePoints, 2);

forceBeta = betaList(idx);
forceGamma = gammaList(idx);

[ay1, ~] = calc_rft_alpha(0, 0, sf);

ay1 = zeros(1,numofForce) + ay1;
 
[ax23, az23] = calc_rft_alpha(forceBeta, forceGamma, sf);
magF1 = -ay1 .* depthList .* (areaList(idx) .* 10 ^ -3);
magF2 = -ax23 .* depthList .* (areaList(idx) .* 10 ^ -3);
% F1 force in N
F1tilde = magF1 .* e1List(:,idx);
% F2 force in N
F2tilde = magF2 .* e2List(:,idx);
F2tilde(3,:) = az23 .* depthList .* (areaList(idx) .* 10^-3);

% scaling factor for angles
phiList = phi(idx);
[f1List,f2List] = normalScale(phiList);

netForce = f1List .* F1tilde + f2List .* F2tilde;
Force = sum(netForce, 2);
end



function [e1List, e2List, vList, vHoriList, v1List, v23List, phi] = calc_velocity(normalList, pointList, wr, vcenter, radius, slipAngle)

% angular velocity radius/s
vcorx = -vcenter * sin(slipAngle);
vcory = vcenter * cos(slipAngle);
vcorz = 0;
vcor = [vcorx; vcory; vcorz];
w = -wr / radius;

numofNormal = size(normalList, 2);
%e2 unit axis
e2List = [normalList(1, :);
    normalList(2, :);
    zeros(1, numofNormal)];
magne2 = sqrt(e2List(1, :) .^2 + e2List(2, :) .^2 + e2List(3, :) .^2);
e2List = [e2List(1, :) ./ magne2;
    e2List(2, :) ./ magne2;
    e2List(3, :) ./ magne2;];

idx2 = (abs(e2List(1, :)) > abs(e2List(2,:)));
temp2 = e2List(1, idx2);
e2List(1, idx2) = e2List(2, idx2);
e2List(2, idx2) = temp2;
% velocity
rList = sqrt(pointList(2, :) .^ 2 + pointList(3, :) .^ 2);   
angleList = atan2(pointList(3, :), pointList(2, :)) + pi/2;

vx = zeros([1, size(angleList, 2)]) + vcor(1);
vy = cos(angleList) .* rList .* w + vcor(2);
vz = sin(angleList) .* rList .* w + vcor(3);

vList = [vx; 
    vy; 
    vz];

% horizontal velocity
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

% calculate phi
phi = calc_Angles(vHoriList, e2List);
phi = wrapToPi(phi);
end



function [betaList, gammaList] = calc_BetaGamma(normalList, e2List, v23List)

% beta
betaList = calc_Angles(normalList, e2List);
idxBeta = normalList(3,:) < 0;
betaList(idxBeta) = - betaList(idxBeta);
betaList = betaList + pi/2;
betaList = wrapToPi(betaList);

% gamma
gammaList = calc_Angles(v23List, e2List);
idxGamma = v23List(3,:) > 0;
gammaList(idxGamma) = - gammaList(idxGamma);
gammaList = wrapToPi(gammaList);

end


% Compute f1, f23 from v1 v23 v
% the sigmoid functions are from [Treers et al., 2021]
function [f1, f23] = normalScale(phi)
global a1 a2 a3 a4 b1 b2 b3 b4
idx1 = (phi > pi/2 & phi <= pi);
idx2 = (phi < 0 & phi >= -pi/2);
idx3 = (phi < -pi/2 & phi >= -pi);

phi(idx1) = pi - phi(idx1);
phi(idx2) = -phi(idx2);
phi(idx3) = phi(idx3) + pi;

v1v = sin(phi);
v23v = cos(phi);

F1 = a1 .* tanh(a2 .* v1v - a3) + a4;
F23 = b1 .* atanh(b2 .* v23v - b3) + b4;
F1max = a1 .* tanh(a2 .* sin(pi/2) - a3) + a4;
F23max = b1 .* atanh(b2 .* cos(0) - b3) + b4;

f1 = F1 ./ F1max;
f23 = F23 ./ F23max;
end


function angles = calc_Angles(v1, v2)
dotprd =v1(1, :) .* v2(1, :) + v1(2, :) .* v2(2, :) + v1(3, :) .* v2(3, :);
timeprd = sqrt(v1(1, :) .^2 + v1(2, :) .^ 2 ...
    + v1(3,: ) .^2) .* sqrt(v2(1, :) .^ 2 ...
    + v2(2,:).^2+v2(3,:).^2);

% not stable, may output imaginary value when angle is around pi. 
angles = acos(dotprd ./ timeprd);
end


% find the local alphax and alphaz with give gamma and beta
% return alpha in N/(cm^3)
function [alphaX, alphaZ] = calc_rft_alpha(beta,gamma, sf)
% using discrete Fourier transform fitting function [Li et al., 2013]
% Fourier coefficients M

global M
beta = wrapToPi(beta);
gamma = wrapToPi(gamma);

idxb1 = (beta >= -pi & beta <= -pi/2);
idxb2 = (beta >= pi/2 & beta <= pi);

beta(idxb1) = beta(idxb1) + pi;
beta(idxb2) = beta(idxb2) - pi;
beta = wrapToPi(beta);

idxg1 = (gamma >= -pi & gamma <= -pi/2);
idxg2 = (gamma >= pi/2 & gamma <= pi);

beta(idxg1) = -beta(idxg1);
beta(idxg2) = -beta(idxg2);

gamma(idxg1) = -pi - gamma(idxg1);
gamma(idxg2) = pi - gamma(idxg2);

gamma = wrapToPi(gamma);

alphaZ = sf .* (M(1) .* cos(0) ...
    + M(2) .* cos(2 .* beta)...
    + M(3) .* sin(2 .* beta + gamma)...
    + M(4) .* sin(gamma)...
    + M(5) .* sin((-2 .* beta) + gamma));

alphaX = sf .* (M(6) .* cos(2 .* beta + gamma)...
    + M(7) .* cos(gamma)...
    + M(8) .* sin(-2 .* beta + gamma)...
    + M(9) .* sin(2 .* beta));

alphaX(idxg1) = -alphaX(idxg1);
alphaX(idxg2) = -alphaX(idxg2);
end


%UNUSED FUNCTIONS

% function [alphaX, alphaZ] = calc_rft_alpha_old(beta,gamma, sf)
% global M
% 
% beta = wrapToPi(beta);
% gamma = wrapToPi(gamma);
% if beta >= -pi && beta <= -pi/2
%     beta = beta + pi;
% elseif beta >= pi/2 && beta <= pi
%     beta = beta - pi;
% end
% beta = wrapToPi(beta);
% gamma = wrapToPi(gamma);
% if gamma >= -pi && gamma <= -pi/2
%     alphaZ = sf*(M(1)*cos(0)+M(2)*cos(2*(-beta))+M(3)*sin(2*(-beta)+(-pi-gamma))...
%         +M(4)*sin((-pi-gamma))+M(5)*sin((-2*(-beta))+(-pi-gamma)));
%     alphaX = -sf*(M(6)*cos(2*(-beta)+(-pi-gamma))+M(7)*cos((-pi-gamma))...
%         +M(8)*sin(-2*(-beta)+(-pi-gamma))+M(9)*sin(2*(-beta)));
% elseif gamma >= pi/2 && gamma <= pi
%     alphaZ = sf*(M(1)*cos(0)+M(2)*cos(2*(-beta))+M(3)*sin(2*(-beta)+(pi-gamma))...
%         +M(4)*sin((pi-gamma))+M(5)*sin((-2*(-beta))+(pi-gamma)));
%     alphaX = -sf*(M(6)*cos(2*(-beta)+(pi-gamma))+M(7)*cos((pi-gamma))...
%         +M(8)*sin(-2*(-beta)+(pi-gamma))+M(9)*sin(2*(-beta)));
% else
%     alphaZ = sf*(M(1)*cos(0)+M(2)*cos(2*beta)+M(3)*sin(2*beta+gamma)...
%         +M(4)*sin(gamma)+M(5)*sin((-2*beta)+gamma));
%     alphaX = sf*(M(6)*cos(2*beta+gamma)+M(7)*cos(gamma)...
%         +M(8)*sin(-2*beta+gamma)+M(9)*sin(2*beta));
% end
% 
% end


% [e1List, e2List] = calc_e1e2(normalList);
% [vList, vHoriList, v1List, v23List] = calc_Vel(pointList, w, radius, e1List, vcor);

% function [e1List, e2List] = calc_e1e2(normalList)
% numofNormal = size(normalList, 2);
% %e2 unit axis
% e2List = [normalList(1, :);
%     normalList(2, :);
%     zeros(1, numofNormal)];
% idxe2 = (e2List(1, :) == 0 & e2List(2, :) == 0);
% e2List(2,:) = 1 .* idxe2 + e2List(2,:) .* (~idxe2);
% 
% magne2 = sqrt(e2List(1, :) .^2 + e2List(2, :) .^2 + e2List(3, :) .^2);
% e2List = [e2List(1, :) ./ magne2;
%     e2List(2, :) ./ magne2;
%     e2List(3, :) ./ magne2;];
% %e1 unit axis
% e1List = [e2List(1,:).*cos(-pi/2)+e2List(2,:).*-sin(-pi/2);
%     e2List(1,:).*sin(-pi/2)+e2List(2,:).*cos(-pi/2);
%     e2List(3,:)];
% 
% end
% function [vList, vHoriList, v1List, v23List] = calc_Vel(pointList, w, radius,e1List,vcor)
% % element radius to the center information (1*n)
% rList = sqrt(pointList(2,:).^2+pointList(3,:).^2);
% 
% % element tangent angle (1*n)
% angleList = atan2(pointList(3,:),pointList(2,:))+pi/2;
% 
% vx = zeros([1,size(angleList,2)])+vcor(1);
% vy = cos(angleList).*rList.*w+vcor(2);
% vz = sin(angleList).*rList.*w+vcor(3);
% 
% vList = [vx;vy;vz];
% 
% vHoriList = [vx;vy;zeros(1,size(angleList,2))];
% 
% idxV = (vHoriList(1,:)==0 & vHoriList(2,:)==0);
% vHoriList(2,:) = 1.*idxV+vHoriList(2,:).*(~idxV);
% v1List = dot(vList,e1List)./1.*e1List;
% v23List = vList-v1List;
% end
