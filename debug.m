wheeldata = matfile('data/smooth_wheel_125.mat');

pointList = wheeldata.Points;

% pointList(1,:) = pointList(1,:) - 30;
depth = 0.04; %m

slipAngle = 30;

[idx, depthList] = run_extractHmap(pointList, slipAngle, depth);

