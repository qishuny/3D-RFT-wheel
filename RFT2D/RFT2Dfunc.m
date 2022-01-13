function forces = RFT2Dfunc(wheeldata, ang_vel, vcenter, sinkage, radius, wheelWidth, scale)

pointList = wheeldata.Points;
areaList = wheeldata.Area;
normalList = wheeldata.Normals;




[vList, e2List] = calc_velocity(normalList, pointList, ang_vel, vcenter, radius);
[betaList, gammaList] = calc_BetaGamma(normalList, vList, e2List);
[netX, netZ, idx, FxList, FzList] = calc_rft_2d(pointList, betaList, gammaList, areaList, sinkage, radius, scale);
Fx = netX * wheelWidth;
Fz = netZ * wheelWidth;
forces = [Fx, Fz];
 
% figure()
% quiver(pointList(1,:), pointList(2,:), vList(1,:), vList(2,:));
% hold on 
% quiver(pointList(1,:), pointList(2,:), e2List(1,:), e2List(2,:));

% figure()
% quiver(pointList(1,:), pointList(2,:), normalList(1,:), normalList(2,:));
% hold on 
% text(pointList(1,:),pointList(2,:),string(betaList(:)*180/pi))

% figure()
% quiver(pointList(1,:), pointList(2,:), vList(1,:), vList(2,:));
% hold on 
% text(pointList(1,:),pointList(2,:),string(gammaList(:)*180/pi))
 

% figure()
% quiver(pointList(1,idx), pointList(2,idx), FxList(1,:), FzList(1,:));
% daspect([1 1 1])

function [vList, e2List] = calc_velocity(normalList, pointList, ang_vel, vcenter, radius)

% omega = [0; ang_vel; 0]
% ri = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)] * [radius;0;0]
% vi = vel + cross(omega, ri)


numofNormal = size(normalList, 2);


% velocity
rList = sqrt(pointList(1, :) .^ 2 + pointList(2, :) .^ 2);   
angleList = atan2(pointList(2, :), pointList(1, :));

vx = sin(angleList) .* rList .* ang_vel + vcenter(1);
vz = -cos(angleList) .* rList .* ang_vel;

vList = [vx;  
    vz];


%e2 unit axis
e2List = [-ones(1, numofNormal);
    zeros(1, numofNormal)];
end

function [betaList, gammaList] = calc_BetaGamma(normalList, vList, e2List)

% beta
betaList = calc_Angles(normalList, e2List);
idxBeta = normalList(2,:) < 0;
betaList = betaList - pi/2;
betaList(idxBeta) = - betaList(idxBeta);

% betaList = wrapToPi(betaList);

% gamma
gammaList = calc_Angles(vList, e2List);
idxGamma = vList(2,:) > 0;
gammaList(idxGamma) = - gammaList(idxGamma);
% gammaList = wrapToPi(gammaList);

end

function angles = calc_Angles(v1, v2)
dotprd =v1(1, :) .* v2(1, :) + v1(2, :) .* v2(2, :);
timeprd = sqrt(v1(1, :) .^2 + v1(2, :) .^ 2)...
    .* sqrt(v2(1, :) .^ 2 + v2(2,:).^2);

% not stable, may output imaginary value when angle is around pi. 
angles = acos(dotprd ./ timeprd);
end

% find the local alphax and alphaz with give gamma and beta
% return alpha in N/(cm^3)
function [alphaX, alphaZ] = calc_alpha(beta, gamma, sf)
% using discrete Fourier transform fitting function [Li et al., 2013]
% Fourier coefficients M
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

% beta = wrapToPi(beta);
% gamma = wrapToPi(gamma);
% 
% idxb1 = (beta >= -pi & beta <= -pi/2);
% idxb2 = (beta >= pi/2 & beta <= pi);
% 
% beta(idxb1) = beta(idxb1) + pi;
% beta(idxb2) = beta(idxb2) - pi;
% beta = wrapToPi(beta);
% 
% idxg1 = (gamma >= -pi & gamma <= -pi/2);
% idxg2 = (gamma >= pi/2 & gamma <= pi);
% 
% beta(idxg1) = -beta(idxg1);
% beta(idxg2) = -beta(idxg2);
% 
% gamma(idxg1) = -pi - gamma(idxg1);
% gamma(idxg2) = pi - gamma(idxg2);
% 
% gamma = wrapToPi(gamma);

alphaZ = sf .* (M(1) .* cos(0) ...
    + M(2) .* cos(2 .* beta)...
    + M(3) .* sin(2 .* beta + gamma)...
    + M(4) .* sin(gamma)...
    + M(5) .* sin((-2 .* beta) + gamma));

alphaX = sf .* (M(6) .* cos(2 .* beta + gamma)...
    + M(7) .* cos(gamma)...
    + M(8) .* cos(-2 .* beta + gamma)...
    + M(9) .* sin(2 .* beta));

% alphaX(idxg1) = -alphaX(idxg1);
% alphaX(idxg2) = -alphaX(idxg2);
end

function [netX, netZ, idx, FxList, FzList] = calc_rft_2d(pointList, betaList, gammaList, areaList, sinkage, radius, scale)
% depth = sand with respect to the center of the wheel mm
depth = -radius + sinkage;

% find points below the surface of the soil
idx = pointList(2,:) < depth;


forceBeta = betaList(idx);
forceGamma = gammaList(idx);
forceDepth = depth - pointList(2,idx);
forceArea = areaList(idx);


[ax, az] = calc_alpha(forceBeta, forceGamma, scale);
FxList = ax .* forceDepth .* forceArea .* 10 ^ -3;
FzList = az .* forceDepth .* forceArea .* 10 ^ -3;

netX = sum(FxList);
netZ = sum(FzList);

% netForce = [FxList; 
%     FzList];
% Force = [netX; netZ];
end

function [sfList] = calc_sf(depth, coeff1, coeff2)
    depth = depth .* 0.001;
    sfList = coeff1 .* (depth .^ coeff2) ./ depth ./ 1000;
end



end
