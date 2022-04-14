function [forces] = RFT3Dfunc(wheeldata, radius, slipAngle, w, vcenter, sinkage, scale, plot)

if nargin < 8 || isempty(plot)
        plot = 0;
end
pointList = wheeldata.Points';
areaList = wheeldata.Area';
normalList = wheeldata.Normals';

depth = -radius + sinkage;
% calculate velocity of each plate element based on body velocities
n_elements = size(pointList, 1);

% angular velocity radius/s
vcorx = -vcenter * sin(slipAngle);
vcory = vcenter * cos(slipAngle);
vcorz = 0;
vcor = [vcorx; vcory; vcorz];

% velocity
rList = sqrt(pointList(:, 2) .^ 2 + pointList(:, 3) .^ 2);   
angleList = atan2(pointList(:, 3), pointList(:, 2)) + pi/2;

vx = zeros(n_elements,1) + vcor(1);
vy = cos(angleList) .* rList .* -w + vcor(2);
vz = sin(angleList) .* rList .* -w + vcor(3);

v_vec = [vx vy vz];
v_norm_vec = v_vec ./ vecnorm(v_vec, 2, 2);
n_norm_vec = normalList ./ vecnorm(normalList, 2, 2);
forces = 0;

% logical conditions that determine whether or not to count force on element
leading_edge = dot(n_norm_vec, v_norm_vec, 2) > 0;      % plate's normal vector has a component in the penetration direction, so it maintains contact forces within sand grains
intruding = pointList(:,3) < depth;                                  % plate is below sand surface
include = leading_edge & intruding;

depthList = depth - pointList(:,3);
% isolate the elements that satisfy this condition for calculation
n_inc = n_norm_vec(include,:);
v_inc = v_norm_vec(include,:);
a_inc = areaList(include,:);
c_inc = pointList(include,:);
d_inc = depthList(include,:);


% generate basis for RFT decomposition
% edge case: if plate is roughly horizontal, use v to define e2 -> stick to 2D RFT model
e2_vec = (abs(n_inc(:,3)) >= 1 - 1e-3) .* ((v_inc+[1e-5 0 0]) .* [1 1 0]) ./ vecnorm((v_inc+[1e-5 0 0]) .* [1 1 0], 2, 2) + ...  % edge case
         (abs(n_inc(:,3)) <  1 - 1e-3) .* (n_inc .* [1 1 0]) ./ vecnorm(n_inc .* [1 1 0], 2, 2);                                 % nominal case
e1_vec = cross(e2_vec, repmat([0,0,1], sum(include), 1));

% decompose velocity vector
v1_vec = dot(v_inc, e1_vec, 2) .* e1_vec;
v23_vec = v_inc - v1_vec;

% calculate 2D RFT parameters
beta_vec =   atan2(  n_inc(:,3), dot(  n_inc, e2_vec, 2)) + pi/2;
gamma_vec = -atan2(v23_vec(:,3), dot(v23_vec, e2_vec, 2));

[aX, aZ] = calc_rft_alpha(beta_vec, gamma_vec, scale);
% 3D alpha component, uses model from 2021 3D RFT
% simplest form: alpha_y = alpha_x(beta=0, gamma=0)
% linear fit form:
% aY = -0.0015 * beta_vec * pi/180 + 0.2035;
aY = calc_rft_alpha(0, 0, scale);
% calculate scaling factors due to decomposition of velocity into 
% tangential (v1) and normal (v23) directions
vt = vecnorm(v1_vec');
ct_fit = [0.440850096954369, 3.62263982880971, 1.60910808139526, 0.832324700111401];
f1 =  ct_fit(1) * ( tanh(ct_fit(2)*vt - ct_fit(3)) +  tanh(ct_fit(3))) / ct_fit(4);

vn = vecnorm(v23_vec');    
cn_fit = [1.99392673405210, 1.61146827229181, 0.973746396532650, 5.81083560370960];
f23 = cn_fit(1) * (atanh(cn_fit(2)*vn - cn_fit(3)) + atanh(cn_fit(3))) / cn_fit(4);
%%%%%%%%%%%% END: MODELS FROM EMPIRICAL DATA %%%%%%%%%%%%

% calculate forces and return vectors to inertial frame
F1 =  -f1' .* aY .* e1_vec .* sign(dot(v_inc, e1_vec, 2));
F2 = -f23' .* aX .* e2_vec;
F3 =  f23' .* aZ .* [0,0,1];
F_i = F1 + F2 + F3;

% superposition: sum force from each plate element to calculate total force & moment
Fi_mat = d_inc .* a_inc .* F_i .* 10^-3;
Fi_mag_mat = vecnorm(Fi_mat');
forces = sum(Fi_mat)';
forces(1) = - forces(1);


if plot == 1
    figure
    
    plot3(pointList(:,1),pointList(:,2),pointList(:,3),'.','Color',[0.6,0.6,0.6],'MarkerSize',1)
    pointList1 = pointList(include,:);
    hold on 
    scatter3(pointList1(:,1),pointList1(:,2),pointList1(:,3), 5, Fi_mag_mat', 'filled','o')
%     title('3D-RFT Forces on Grousered Wheel');
    daspect([1 1 1])
    cb = colorbar('southoutside'); 
    cb.FontSize = 16;
    view(-55,25)
    axis off
    
    figure
    hold on
    pointList1 = pointList(include,:);
    X = pointList1(:,1);
    Y = pointList1(:,2);
    Z = pointList1(:,3);
    U = Fi_mat(:,1);
    V = Fi_mat(:,2);
    W = Fi_mat(:,3);
    
    q = quiver3(X, Y, Z, U, V, W, 2);
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
    
    daspect([1 1 1])
    view(-55,25)
    axis off
    
%     figure
%     quiver3(c_inc(:,1),c_inc(:,2),c_inc(:,3),n_inc(:,1),n_inc(:,2),n_inc(:,3),1,'Color', [0,0.2,0.8]);
% 
%     figure
%     quiver3(c_inc(:,1),c_inc(:,2),c_inc(:,3),v_inc(:,1),v_inc(:,2),v_inc(:,3),1,'Color', [0,0.2,0.8]);

end

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

end