

function [idxOut, depthList, pile, under] = run_extractHmap(pointList, slipAngle, sink, plotToggle)

pointList(1,:) = pointList(1,:) - 0.5*(max(pointList(1,:))-min(pointList(1,:)));
pointListOriginal = pointList;
num = size(pointList, 2);

wheelDiameter = 0.125; %m
wheelWidth = 0.06; %m
sinkage = 0.04; %m

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
exclude = z <= - 100;
f1 = fit([x y],z,'poly55', 'Exclude', exclude);

for i = 1:size(Xtrimed,1)
    for j = 1:size(Xtrimed,2)
        
        SandHmapnew(i, j) = f1(Xtrimed(i, j), Ytrimed(i, j));       
    end
end



pointListOriginal = rotateZ(pointListOriginal, -slipAngle * pi/180);
spz = interp2(double(Xtrimed), double(Ytrimed), double(SandHmapnew), pointListOriginal(1,:)', pointListOriginal(2,:)');
spz = spz';

depth =  - wheelDiameter/ 2 + sink * 1000;
depth
pile = pointListOriginal(3,:) < spz & pointListOriginal(3,:) > depth;
under = pointListOriginal(3,:) < spz & pointListOriginal(3,:) <= depth;
idxOut = pile | under;

depthList = zeros(1, num);
depthList(pile) = abs(spz(pile) - pointList(3,pile)) * 0.1693;

depthListPile = depthList(pile);
pointListPile = pointList(:,pile);

depthList(under) = min(abs(depth - pointList(3,under)), abs(spz(under) - pointList(3,under)));
depthListUnder = depthList(under);
pointListUnder = pointList(:,under);

depthList(~idxOut) = 0;




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
    
    figure
    plot3(pointList(1,pile), pointList(2,pile), pointList(3,pile),'.','Color',[0.1,0.1,0.1],'MarkerSize',3);
    hold on
    plot3(pointList(1,under), pointList(2,under), pointList(3,under),'.','Color',[0.2,0.2,0.2],'MarkerSize',3);
    plot3(pointList(1,(~under & ~pile)), pointList(2,(~under & ~pile)), pointList(3,(~under & ~pile)),'.','Color',[0.9,0.9,0.9],'MarkerSize',1);
    s = surf(Xtrimed, Ytrimed, SandHmapnew, 'FaceAlpha', 0.4);
    s.EdgeColor = [0.9 0.9 0.9];
    s.FaceColor = [1 1 0];
    axis equal 
    daspect([1 1 1])
    view(-105,25)
    figure
    plot3(pointList(1,pile), pointList(2,pile), pointList(3,pile),'.','Color','r','MarkerSize',7);
    hold on

    
    plot3(pointList(1,under), pointList(2,under), pointList(3,under),'.','Color','b','MarkerSize',7);
    
    plot3(pointList(1,(~under & ~pile)), pointList(2,(~under & ~pile)), pointList(3,(~under & ~pile)),'.','Color',[0.8,0.8,0.8],'MarkerSize',1);
    s = surf(Xtrimed, Ytrimed, SandHmapnew, 'FaceAlpha', 0.5);
    s.EdgeColor = 'none';
    axis equal 
    legend('pile-up','undisturbed')
    view(-105,25)

   
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