% wheeldata = matfile('data/smooth_wheel_125.mat');
% 
% pointList = wheeldata.Points;
% 
% % pointList(1,:) = pointList(1,:) - 30;
% depth = 0.04; %m
% 
% slipAngle = 30;
% 
% [idx, depthList] = run_extractHmap(pointList, slipAngle, depth);
v1 = [0.011764507331288; 0.011764507331288; 0.676209435691110]
v2 = [-0.707106781186548; -0.707106781186548; 0]


% angles = calc_Angles(v23List(:,522), e2List(:,522))
angles = calc_Angles(v23List, e2List);
function angles = calc_Angles(v1, v2)
dotprd =v1(1, :) .* v2(1, :) + v1(2, :) .* v2(2, :) + v1(3, :) .* v2(3, :);
timeprd = sqrt(v1(1, :) .^2 + v1(2, :) .^ 2 ...
    + v1(3,: ) .^2) .* sqrt(v2(1, :) .^ 2 ...
    + v2(2,:).^2+v2(3,:).^2);

% not stable, may output imaginary value when angle is around pi. 
angles = acos(dotprd ./ timeprd);
angles = real(angles);
end