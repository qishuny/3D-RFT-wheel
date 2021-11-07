
wheeldata = matfile('data/smooth_wheel_125.mat');
pointList = wheeldata.Points;

pointList(1,:) = pointList(1,:) - 0.5*(max(pointList(1,:))-min(pointList(1,:)));
pointListOriginal = pointList;
% slipAngle in degree
slipAngle = 45;

wheelDiameter = 0.125; %m
wheelWidth = 0.06; %m
depth = 0.04; %m

%grid per m
n = 200;

gridsize = 1/n;

[pointList] = rotateZ(pointList, -slipAngle * pi / 180);


% 1 for plot all values
plotToggle = 0;

tic
[sandHmap, wheelPos] = extractHmap((90 - slipAngle), wheelDiameter, wheelWidth, depth, n);
sandHmap = sandHmap';
toc

%% change unit to mm
sandHmap = sandHmap * 1000;
wheelPos = wheelPos * 1000;
wheelDiameter = wheelDiameter * 1000;
wheelWidth = wheelWidth * 1000;
depth = depth * 1000; 
gap = gridsize * 1000;
%% plot simluation result: sand height map & wheel position
close all

[X, Y] = meshgrid(-400:gap:400, -400:gap:400);

surf(X, Y, sandHmap, 'FaceAlpha', 0.5)
%     plot simluation result: sand height map & wheel position
figure
s = surf(X, Y, sandHmap);
s.EdgeColor = 'none';
hold on
scatter3(wheelPos(1), wheelPos(2), wheelPos(3),'r')
axis on
xlabel('x')
ylabel('y')
axis equal
%% line up the wheel and the sand height map

lim = sqrt((wheelDiameter / 2) ^ 2 + (wheelWidth / 2) ^ 2) + 10;
idxX = X(1, :) >(wheelPos(1) - lim) & X(1, :) < (wheelPos(1) + lim);
idxY = Y(:, 1) >(wheelPos(2) - lim) & Y(:, 1) < (wheelPos(2) + lim);

Y = Y - wheelPos(2);
sandHmap = sandHmap - wheelPos(3);

wheelPos(3) = 0;

Xtrimed = X(idxY, idxX);
Ytrimed = Y(idxY, idxX);
SandHmapOriginal = sandHmap(idxY, idxX);




%% 


XY = [Xtrimed(:) Ytrimed(:)];  
theta = slipAngle;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
rotXY = XY * R';
Xqr = reshape(rotXY(:,1), size(Xtrimed,1), []);
Yqr = reshape(rotXY(:,2), size(Ytrimed,1), []);






idxWheelX = Xqr(:, :) <=  wheelWidth / 2 + 2 & Xqr(:, :) >= - wheelWidth/ 2 - 2;
idxWheelY = Yqr(:, :) <= wheelDiameter/ 2  + 2& Yqr(:, :) >= - wheelDiameter/ 2 -2;

idx = idxWheelX & idxWheelY;
SandHmapnew = SandHmapOriginal;

% SandHmapnew(idx) = max(SandHmapnew(idx), -depth);
SandHmapnew(idx) = -100;
% SandHmapnew(idx) = -depth;


x = reshape(Xtrimed,[],1);
y = reshape(Ytrimed,[],1);
z = reshape(SandHmapnew,[],1);
f1 = fit([x y],z,'lowess');

for i = 1:size(Xtrimed,1)
    for j = 1:size(Xtrimed,2)
        
        SandHmapnew(i, j) = f1(Xtrimed(i, j), Ytrimed(i, j));       
    end
end

figure
plot(f1, [x y], z);

figure
s = surf(Xtrimed, Ytrimed, SandHmapnew, 'FaceAlpha', 0.5);
title("sand height map of slip angle", slipAngle)
pointListOriginal = rotateZ(pointListOriginal, -slipAngle * pi/180);
spz = interp2(Xtrimed, Ytrimed, SandHmapnew, pointListOriginal(1,:)', pointListOriginal(2,:)');
spz = spz';

under = pointListOriginal(3,:) < spz;
depthList = abs(spz(under) - pointListOriginal(3,under));

figure
plot3(pointListOriginal(1,under), pointListOriginal(2,under), pointListOriginal(3,under),'.','Color','r','MarkerSize',1);
hold on
plot3(pointListOriginal(1,~under), pointListOriginal(2,~under), pointListOriginal(3,~under),'.','Color',[0.8,0.8,0.8],'MarkerSize',1);
s = surf(Xtrimed, Ytrimed, SandHmapnew, 'FaceAlpha', 0.5);
s.EdgeColor = 'none';
axis equal








if plotToggle == 1
    surf(X, Y, sandHmap, 'FaceAlpha', 0.5)
    % plot simluation result: sand height map & wheel position
    figure
    s = surf(X, Y, sandHmap);
    s.EdgeColor = 'none';
    hold on
    scatter3(wheelPos(1), wheelPos(2), wheelPos(3),'r')
    axis on
    xlabel('x')
    ylabel('y')
    axis equal

    
    % plot orignial sand height map
    figure

    plot3(pointList(1,:), pointList(2,:), pointList(3,:), '.', ...
        'Color', [0.0,0.0,0.0], 'MarkerSize', 1)
    hold on
    surf(Xtrimed, Ytrimed, SandHmapnew, 'FaceAlpha', 0.5)

    axis on
    xlabel('x')
    ylabel('y')
    axis equal
    
    % plot inflated wheel and modified sand height map
    figure
    plot3(Points_inflated(1,under), Points_inflated(2,under), Points_inflated(3,under),'.','Color','r','MarkerSize',1);
    hold on
    plot3(Points_inflated(1,~under), Points_inflated(2,~under), Points_inflated(3,~under),'.','Color',[0.8,0.8,0.8],'MarkerSize',1);
    s = surf(Xtrimed, Ytrimed, SandHmapnew, 'FaceAlpha', 0.5);
    s.EdgeColor = 'none';
    axis equal
    
    % check depth of the wheel
%     forcePoints = Points_inflated(:, under);
%     figure
%     testG = 100;
%     for i = 1:testG:size(forcePoints,2)
%         plot3(forcePoints(1,i), forcePoints(2,i), forcePoints(3,i),'.','Color','r','MarkerSize',1);
%         hold on
%         text(forcePoints(1,i),forcePoints(2,i),forcePoints(3,i),string(depthList(i)));
% 
%     end
%     s = surf(Xtrimed, Ytrimed, SandHmapnew, 'FaceAlpha', 0.5);
    
    % plot wheel points above & below the sand
   
end


function [Points] = rotateZ(Points,thetaz)

Rz = [cos(thetaz), -sin(thetaz), 0;
    sin(thetaz), cos(thetaz), 0;
    0, 0, 1;];
Points = Rz * Points;

end

function [Points_inflated] = inflate(Points, wheelDiameter, wheelWidth, gap)
scaleY = (wheelDiameter / 2 + gap) / (wheelDiameter /2);
scaleX = (wheelWidth / 2 + gap) / (wheelWidth /2);

Points_inflated(1,:) = Points(1,:) .* scaleX;
Points_inflated(2,:) = Points(2,:) .* scaleY;
Points_inflated(3,:) = Points(3,:);
end 