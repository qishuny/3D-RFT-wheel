% 3D RFT Modeling on Wheels
% Qishun Yu Catherine Pavlov
% 02/21/2021
function [forces] = RFT3Dfunc(wheeldata, radius, slipAngle, wr, vcenter, sinkage, scale)
% load wheel point data
pointList = wheeldata.Points;
areaList = wheeldata.Area;
normalList = wheeldata.Normals;

if wr == 0
    wr = 0.00001;
end


%% Single RFT calc

% Geometry & Velocity Calc

[e1List, e2List, vList, vHoriList, v1List, v23List, phi] = calc_velocity(normalList, pointList, wr, vcenter, radius, slipAngle);
% Force Calc
[betaList, gammaList] = calc_BetaGamma(normalList, e2List, v23List);
[Fx, Fy, Fz, idx, netForce] = calc_3D_rft(pointList, normalList, vList, betaList, gammaList, e1List, e2List, areaList, phi, sinkage, radius, scale);

forces = [Fx; -Fy; Fz];

plotForce = 0;
plotVelocity = 0 ;
%% plot force
if plotForce == 1
    figure

    plot3(pointList(1,:),pointList(2,:),pointList(3,:),'.','Color',[0.6,0.6,0.6],'MarkerSize',1)

    hold on
    pointList1 = pointList(:,idx);
    X = pointList1(1,:);
    Y = pointList1(2,:);
    Z = pointList1(3,:);
    U = netForce(1,:);
    V = netForce(2,:);
    W = netForce(3,:);
    
    q = quiver3(X, Y, Z, U, V, W);
    mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
            reshape(q.WData, numel(q.UData), [])).^2, 2));
    

    currentColormap = colormap(gca);
    [~, ~, ind] = histcounts(mags, size(currentColormap, 1));

    cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
    cmap(:,:,4) = 255;
    cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

    set(q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'
    set(q.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');
    
    title('3D-RFT Forces on Grousered Wheel');
    daspect([1 1 1])
    view(-55,15)
    axis off
end

%% plot velocity
gapSize = 20;
if plotVelocity == 1  
    %Plot selected velocity and v1,v23
    figure
    

    quiver3(pointList(1,idx),pointList(2,idx),pointList(3,idx),vList(1,idx),vList(2,idx),vList(3,idx),1,'Color', [0,0.2,0.8]);
    hold on
    quiver3(pointList(1,idx),pointList(2,idx),pointList(3,idx),v1List(1,idx),v1List(2,idx),v1List(3,idx),1,'r');
    quiver3(pointList(1,idx),pointList(2,idx),pointList(3,idx),v23List(1,idx),v23List(2,idx),v23List(3,idx),1,'g');
    legend('velocity','v1','v23')
    daspect([1 1 1])
    

    % Plot v23 and e2
%     figure
%     for k =1:gapSize:size(pointList,2)
%         quiver3(pointList(1,k),pointList(2,k),pointList(3,k),e2List(1,k),e2List(2,k),e2List(3,k),'g');
%         hold on
%         quiver3(pointList(1,k),pointList(2,k),pointList(3,k),v23List(1,k),v23List(2,k),v23List(3,k),'r');
%         legend('e2','v23')
%     %     text(pointList(1,k),pointList(2,k),pointList(3,k),string(gammaList(k)*180/pi))
%         daspect([1 1 1])
%     end

    %Plot v1 and e1
%     figure
%     
%     idxS = pointList(3,:) < -62.49;
%     plot3(pointList(1,idxS),pointList(2,idxS),pointList(3,idxS),'ok','MarkerFaceColor',[0,0.5,0.5])
%     for k =1:gapSize:size(pointList,2)
%         quiver3(pointList(1,k),pointList(2,k),pointList(3,k),e1List(1,k),e1List(2,k),e1List(3,k),'g');
%         hold on
%         quiver3(pointList(1,k),pointList(2,k),pointList(3,k),v1List(1,k),v1List(2,k),v1List(3,k),'r');
%         legend('e1','v1')
%         daspect([1 1 1])
%     end
end



function [Fx, Fy, Fz, idx, netForce] = calc_3D_rft(pointList, normalList, vList, betaList, gammaList, e1List, e2List, areaList, phi, sinkage, radius, sf)
    
% depth = sand with respect to the center of the wheel mm
depth = -radius + sinkage;

% find points below the surface of the soil
idx1 = pointList(3,:) < depth;
idx2 = dot(normalList,vList) >= 0;
idx = idx1 & idx2;

forcePoints = pointList(:, idx);
numofForce = size(forcePoints, 2);

forceBeta = betaList(idx);
forceGamma = gammaList(idx);

[ay1, ~] = calc_rft_alpha(0, 0, sf);

ay1 = zeros(1,numofForce) + ay1;
 
[ax23, az23] = calc_rft_alpha(forceBeta, forceGamma, sf);
magF1 = ay1 .* abs(depth - pointList(3,idx)) .* (areaList(idx) .* 10 ^ -3);
magF2 = ax23 .* abs(depth - pointList(3,idx)) .* (areaList(idx) .* 10 ^ -3);

% F1 force in N
F1tilde = magF1 .* e1List(:,idx);
% F2 force in N
F2tilde = magF2 .* e2List(:,idx);
F2tilde(3,:) = az23 .* abs(depth - pointList(3,idx)) .* (areaList(idx) .* 10^-3);


% scaling factor for angles
phiList = phi(idx);
[f1List,f2List] = normalScale(phiList);

netForce = f1List .* F1tilde + f2List .* F2tilde;
Force = sum(netForce, 2);
Fx = Force(1);
Fy = Force(2);
Fz = Force(3);
end



function [e1List, e2List, vList, vHoriList, v1List, v23List, phi] = calc_velocity(normalList, pointList, w, vcenter, radius, slipAngle)

% angular velocity radius/s
vcorx = -vcenter * sin(slipAngle);
vcory = vcenter * cos(slipAngle);
vcorz = 0;
vcor = [vcorx; vcory; vcorz];


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
vy = cos(angleList) .* rList .* -w + vcor(2);
vz = sin(angleList) .* rList .* -w + vcor(3);

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
a1 = 0.44;
a2 = 3.62;
a3 = 1.61;
a4 = 0.41;
b1 = 1.99;
b2 = 1.61;
b3 = 0.97;
b4 = 4.31;

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

function [sfList] = calc_sf(depth, coeff1, coeff2)
    depth = depth .* 0.001;
    sfList = coeff1 .* (depth .^ coeff2) ./ depth ./ 1000;
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
function [alphaX, alphaZ] = calc_rft_alpha(beta, gamma, sf)
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
end
