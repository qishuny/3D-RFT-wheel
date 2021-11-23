

function [idxOut, depthListOut, pile, under] = DEMdepthList(pointList, slipAngle, depth, sf1, sf2)

slipAngle_degree = slipAngle*180/pi;
num = size(pointList, 2);
pointList(1,:) = pointList(1,:) - 0.5*(max(pointList(1,:))-min(pointList(1,:)));
filename = strcat('output/sand', num2str(slipAngle_degree), '.mat');

sandData = matfile(filename);
Xtrimed = sandData.Xtrimed;

Ytrimed = sandData.Ytrimed;
SandHmapnew = sandData.SandHmapnew;

pointList = rotateZ(pointList, -slipAngle);
spz = interp2(Xtrimed, Ytrimed, SandHmapnew, pointList(1,:)', pointList(2,:)');
spz = spz';



pile = pointList(3,:) < spz & pointList(3,:) > depth;
under = pointList(3,:) < spz & pointList(3,:) <= depth;

idxOut = pile | under;


depthList = zeros(1, num);
depthList(pile) = abs(spz(pile) - pointList(3,pile)) .* sf1;

depthListPile = depthList(pile);
pointListPile = pointList(:,pile);

depthList(under) = min(abs(depth - pointList(3,under)), abs(spz(under) - pointList(3,under))) .* sf2;
% depthList(under) = abs(depth - pointList(3,under));
depthListUnder = depthList(under);
pointListUnder = pointList(:,under);

depthListOut = depthList(idxOut);
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
    
    figure
    plot3(pointList(1,pile), pointList(2,pile), pointList(3,pile),'.','Color','r','MarkerSize',1);
    hold on
%     for k =1:100:size(pointListPile,2)
%         hold on
%         text(pointListPile(1,k),pointListPile(2,k),pointListPile(3,k),string(depthListPile(k)));
%     end
    
    plot3(pointList(1,under), pointList(2,under), pointList(3,under),'.','Color','y','MarkerSize',1);
    
%     for k =1:100:size(depthListUnder,2)
%         hold on
%         text(pointListUnder(1,k),pointListUnder(2,k),pointListUnder(3,k),string(depthListUnder(k)));
%     end
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