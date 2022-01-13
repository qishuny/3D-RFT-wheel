wheeldata = matfile('data/smooth_wheel_125.mat');
% wheeldata = matfile('data/grousered_wheel_125.mat');
pointList = wheeldata.Points;

% pointList(1,:) = pointList(1,:) - 30;
% depth = 0.04; %m


slipAngle = 0;
% SET velocity of the center of rotation of the body mm/s
vcenter = 10;
% SET wheel rotational speed mm/s
wr = 10;
% SET sinkage mm


sf1 = 0.175;
sf2 = 0.5;
flist = zeros(1000, 3);
sinkage = 35.68;
[Fx, Fy, Fz] = RFT3Dfunc(wheeldata, slipAngle, wr, 20, sinkage, sf1, sf2)
% for i = 1:1000
%     sinkage = 10 + i/100*4;
%     [flist(i,1), flist(i,2), flist(i,3)] = RFT3Dfunc(wheeldata, slipAngle, wr, vcenter, sinkage, sf1, sf2);
% end

% [idxOut, depthList, pile, under] = run_extractHmapFitTest(pointList, slipAngle, depth);
% size(depthList)
% [idx, depthList] = run_extractHmap(pointList, slipAngle, depth);
% v1 = [0.011764507331288; 0.011764507331288; 0.676209435691110]
% v2 = [-0.707106781186548; -0.707106781186548; 0]
% 
% 
% % angles = calc_Angles(v23List(:,522), e2List(:,522))
% angles = calc_Angles(v23List, e2List);
% function angles = calc_Angles(v1, v2)
% dotprd =v1(1, :) .* v2(1, :) + v1(2, :) .* v2(2, :) + v1(3, :) .* v2(3, :);
% timeprd = sqrt(v1(1, :) .^2 + v1(2, :) .^ 2 ...
%     + v1(3,: ) .^2) .* sqrt(v2(1, :) .^ 2 ...
%     + v2(2,:).^2+v2(3,:).^2);
% 
% % not stable, may output imaginary value when angle is around pi. 
% angles = acos(dotprd ./ timeprd);
% angles = real(angles);
% end
% coeff1 = 1.915315192989908e+03;
% coeff2 = 1.364833809455482;
% 
% depth = linspace(0, 80, 1000);
% sfList = calc_sf(depth, coeff1, coeff2);
% sf = mean(sfList)
% figure()
% plot(depth, sfList)
% function [sfList] = calc_sf(depth, coeff1, coeff2)
%     depth = depth .* 0.001;
%     sfList = coeff1 .* (depth .^ coeff2) ./ depth ./ 1000;
% end