wheeldata = matfile('smooth_wheel_125.mat');
% wheeldata = matfile('grousered_wheel_125.mat');
pointList = wheeldata.Points';
addpath('../BK_terra') 
radius = 62.5;
vcenter = 10;
scale = 0.6;
sinkage = 40;
wr = 10;
save_extractHmap(pointList', 30, sinkage / 1000, wr, 1)

function save_extractHmap(pointList, slipAngle, sink, wr, plotToggle)

pointList(1,:) = pointList(1,:) - 0.5*(max(pointList(1,:))-min(pointList(1,:)));
pointListOriginal = pointList;
num = size(pointList, 2);

wheelDiameter = 0.125; %m
wheelWidth = 0.06; %m
sinkage = sink; %m

%grid per m
n = 200;

gridsize = 1/n;

[pointList] = rotateZ(pointList, -slipAngle * pi / 180);


tic
[sandHmap, wheelPos] = extractHmap((90 - slipAngle), wheelDiameter, wheelWidth, sinkage, n, plotToggle);
sandHmap = sandHmap';
toc

%% change unit to mm
sandHmap = sandHmap * 1000;
wheelPos = wheelPos * 1000;
wheelDiameter = wheelDiameter * 1000;
wheelWidth = wheelWidth * 1000;
sinkage = sinkage * 1000; 
gap = gridsize * 1000;
%% plot simluation result: sand height map & wheel position

[X, Y] = meshgrid(-400:gap:400, -400:gap:400);


% %% line up the wheel and the sand height map

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






idxWheelX = Xqr(:, :) <=  wheelWidth / 2 + 5 & Xqr(:, :) >= - wheelWidth/ 2 - 5;
idxWheelY = Yqr(:, :) <= wheelDiameter/ 2  + 5& Yqr(:, :) >= - wheelDiameter/ 2 -5;

idx = idxWheelX & idxWheelY;
SandHmapnew = SandHmapOriginal;

% SandHmapnew(idx) = max(SandHmapnew(idx), -depth);
SandHmapnew(idx) = -100;
% SandHmapnew(idx) = -depth;


x = reshape(Xtrimed,[],1);
y = reshape(Ytrimed,[],1);
z = reshape(SandHmapnew,[],1);
exclude = z <= - 100;
f1 = fit([x y],z,'poly55', 'Exclude', exclude);

for i = 1:size(Xtrimed,1)
    for j = 1:size(Xtrimed,2)
        
        SandHmapnew(i, j) = f1(Xtrimed(i, j), Ytrimed(i, j));       
    end
end

fname = sprintf('%dsand%d.mat',slipAngle, wr)
save(fname, 'SandHmapnew', 'Xtrimed', 'Ytrimed');
if plotToggle == 1

    figure
    s1 = surf(Xtrimed, Ytrimed, SandHmapOriginal, 'FaceAlpha', 0.5);
    hold on
    scatter3(wheelPos(1), wheelPos(2), wheelPos(3),'r')
    
    axis equal
    figure
    plot(f1, [x y], z, 'Exclude', exclude);

    figure
    s = surf(Xtrimed, Ytrimed, SandHmapnew, 'FaceAlpha', 0.5);

   
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

end