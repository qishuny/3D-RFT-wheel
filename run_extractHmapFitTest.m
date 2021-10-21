

function [idxOut, depthList, pile, under] = run_extractHmapFitTest(pointList, slipAngle, depth)

num = size(pointList, 2);
pointList(1,:) = pointList(1,:) - 0.5*(max(pointList(1,:))-min(pointList(1,:)));
filename = strcat('output/sand', num2str(90), '.mat');

sandData = matfile(filename);
Xtrimed = sandData.Xtrimed;

Ytrimed = sandData.Ytrimed;
SandHmapnew = sandData.SandHmapnew;

pointList = rotateZ(pointList, -slipAngle * pi/180);
spz = interp2(Xtrimed, Ytrimed, SandHmapnew, pointList(1,:)', pointList(2,:)');
spz = spz';



pile = pointList(3,:) < spz & pointList(3,:) > depth;
under = pointList(3,:) < spz & pointList(3,:) <= depth;

idxOut = pile | under;


depthList = zeros(1, num);
depthList(pile) = abs(spz(pile) - pointList(3,pile)) * 1;
depthList(pile) = 0;
depthList(under) = min(abs(depth - pointList(3,under)), abs(spz(under) - pointList(3,under)));

plotToggle = 0;


if plotToggle == 1
% 
%     figure
%     s = surf(Xtrimed, Ytrimed, SandHmapnew, 'FaceAlpha', 0.5);
    
    figure
    plot3(pointList(1,pile), pointList(2,pile), pointList(3,pile),'.','Color','r','MarkerSize',1);
    hold on
    plot3(pointList(1,under), pointList(2,under), pointList(3,under),'.','Color','y','MarkerSize',1);
    plot3(pointList(1,(~under & ~pile)), pointList(2,(~under & ~pile)), pointList(3,(~under & ~pile)),'.','Color',[0.8,0.8,0.8],'MarkerSize',1);
    s = surf(Xtrimed, Ytrimed, SandHmapnew, 'FaceAlpha', 0.5);
    s.EdgeColor = 'none';
    axis equal 
end


function [Points] = rotateZ(Points,thetaz)

Rz = [cos(thetaz), -sin(thetaz), 0;
    sin(thetaz), cos(thetaz), 0;
    0, 0, 1;];
Points = Rz * Points;

end

end