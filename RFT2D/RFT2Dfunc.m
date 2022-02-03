function forces = RFT2Dfunc(wheeldata, ang_vel, vcenter, sinkage, radius, wheelWidth, scale)
pointList = wheeldata.Points;
areaList = wheeldata.Area;
normalList = wheeldata.Normals;

[vList, e2List] = calc_velocity(normalList, pointList, ang_vel, vcenter, radius);
[betaList, gammaList] = calc_BetaGamma(normalList, vList, e2List);
[netX, netZ, idx, FxList, FzList] = calc_rft_2d(pointList, vList, betaList, gammaList, areaList, sinkage, radius, scale);
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

% gamma
gammaList = calc_Angles(vList, e2List);
idxGamma = vList(2,:) > 0;
gammaList(idxGamma) = - gammaList(idxGamma);
end

function angles = calc_Angles(v1, v2)
dotprd =v1(1, :) .* v2(1, :) + v1(2, :) .* v2(2, :);
timeprd = sqrt(v1(1, :) .^2 + v1(2, :) .^ 2)...
    .* sqrt(v2(1, :) .^ 2 + v2(2,:).^2);
angles = acos(dotprd ./ timeprd);
end

% find the local alphax and alphaz with give gamma and beta
% return alpha in N/(cm^3)
function [alphaX, alphaZ] = calc_alpha(beta, gamma, sf)
% using discrete Fourier transform fitting function [Li et al., 2013]
% Fourier coefficients M
M = [0.206, 0.169, 0.212, 0.358, 0.055, -0.124, 0.253, 0.007, 0.088];

alphaZ = sf .* (M(1) .* cos(0) ...
    + M(2) .* cos(2 .* beta)...
    + M(3) .* sin(2 .* beta + gamma)...
    + M(4) .* sin(gamma)...
    + M(5) .* sin((-2 .* beta) + gamma));

alphaX = sf .* (M(6) .* cos(2 .* beta + gamma)...
    + M(7) .* cos(gamma)...
    + M(8) .* cos(-2 .* beta + gamma)...
    + M(9) .* sin(2 .* beta));
end

function [netX, netZ, idx, FxList, FzList] = calc_rft_2d(pointList, vList, betaList, gammaList, areaList, sinkage, radius, scale)
% depth = sand with respect to the center of the wheel mm
depth = -radius + sinkage;

% find points below the surface of the soil
idx1 = pointList(2,:) < depth;
idx2 = dot(pointList,vList) >= -1e-5;
idx = idx1 & idx2;

forceBeta = betaList(idx);
forceGamma = gammaList(idx);
forceDepth = depth - pointList(2,idx);
forceArea = areaList(idx);

[ax, az] = calc_alpha(forceBeta, forceGamma, scale);
FxList = ax .* forceDepth .* forceArea .* 10 ^ -3;
FzList = az .* forceDepth .* forceArea .* 10 ^ -3;

netX = sum(FxList);
netZ = sum(FzList);
end

end
