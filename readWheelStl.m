clear all
close all
clc

%% grousered wheel
TRG = stlread('grousered_wheel2mm.STL');
trisurf(TRG);
daspect([1 1 1]);

Points_grousered = incenter(TRG)';
Normals_grousered = faceNormal(TRG)';
Area_grousered = generateArea(TRG.Points', TRG.ConnectivityList');


[Points_grousered, Area_grousered, Normals_grousered] = trimGrousered(Points_grousered, Area_grousered, Normals_grousered);
[Points_grousered, Normals_grousered] = rotateX(Points_grousered, Normals_grousered, pi/2);
[Points_grousered, Normals_grousered] = rotateZ(Points_grousered, Normals_grousered, pi/2);
[Points_grousered, Normals_grousered] = rotateZ(Points_grousered, Normals_grousered, pi);
Points_grousered(1,:) = Points_grousered(1,:) + 60;

magneN = sqrt(Normals_grousered(1,:).^2+Normals_grousered(2,:).^2+Normals_grousered(3,:).^2);
Normals_grousered = [Normals_grousered(1,:)./magneN;
    Normals_grousered(2,:)./magneN;
    Normals_grousered(3,:)./magneN;];

figure
quiver3(Points_grousered(1,:),Points_grousered(2,:),Points_grousered(3,:), ...
     Normals_grousered(1,:),Normals_grousered(2,:),Normals_grousered(3,:),1,'color','b');
daspect([1 1 1]);

Points = Points_grousered;
Normals = Normals_grousered;
Area = Area_grousered;
save('grousered_wheel_125.mat','Points','Normals','Area')



%% smooth wheel
TRS = stlread('smooth_wheel2mm.STL');
figure
trisurf(TRS);
daspect([1 1 1]);

Points_smooth = incenter(TRS)';
Normals_smooth = faceNormal(TRS)';  
Area_smooth = generateArea(TRS.Points',TRS.ConnectivityList'); 
[Points_smooth, Area_smooth, Normals_smooth] = trimSmooth(Points_smooth, Area_smooth, Normals_smooth);
% Rotate around x axis
[Points_smooth, Normals_smooth] = rotateX(Points_smooth, Normals_smooth, pi/2);

% Rotate around z axis
[Points_smooth, Normals_smooth] = rotateZ(Points_smooth, Normals_smooth, pi/2);


magneN = sqrt(Normals_smooth(1,:).^2+Normals_smooth(2,:).^2+Normals_smooth(3,:).^2);
Normals_smooth = [Normals_smooth(1,:)./magneN;
    Normals_smooth(2,:)./magneN;
    Normals_smooth(3,:)./magneN;];

figure
quiver3(Points_smooth(1,:),Points_smooth(2,:),Points_smooth(3,:), ...
     Normals_smooth(1,:),Normals_smooth(2,:),Normals_smooth(3,:),1,'color','b');
daspect([1 1 1]);

Points = Points_smooth;
Normals = Normals_smooth;
Area = Area_smooth;
save('smooth_wheel_125.mat','Points','Normals','Area')




function areaarray = generateArea(Points,List)
    areaarray = [];   
    for i = 1:size(List,2)  
        temp = [];
        for j = 1:3            
            idx = List(j,i);     
            temp = [temp [Points(1,idx);Points(2,idx);Points(3,idx)]];

        end
        x1 = temp(1,1);
        x2 = temp(1,2);
        x3 = temp(1,3);
        y1 = temp(2,1);
        y2 = temp(2,2);
        y3 = temp(2,3);
        z1 = temp(3,1);
        z2 = temp(3,2);
        z3 = temp(3,3);
        %r1-r2
        a = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
        %r2-r3
        b = sqrt((x3-x2)^2+(y3-y2)^2+(z3-z2)^2);
        %r1-r3
        c = sqrt((x1-x3)^2+(y1-y3)^2+(z1-z3)^2);
        s = (a+b+c)/2;
        AreaTemp = sqrt(s*(s-a)*(s-b)*(s-c)); 
        areaarray = [areaarray AreaTemp];
    end
end


function [Points, Area, Normals] = trimSmooth(Points, Area, Normals)
% center the wheel to [0,0]
Points(1,:) = Points(1,:)-62.5;
Points(2,:) = Points(2,:)-62.5;

idx1 = Points(3,:)<59.9;

Points = Points(:,idx1);
Area = Area(idx1);
Normals = Normals(:,idx1);

size(Points,2)

radius = (Points(1,:).^2+Points(2,:).^2);
idx2 =(radius<(61^2)) & (Points(3,:)>0);

Points = Points(:,~idx2);
Area = Area(~idx2);
Normals = Normals(:,~idx2);

size(Points,2)
end

function [Points, Area, Normals] = trimGrousered(Points, Area, Normals)
% center the wheel to [0,0]
Points(1,:) = Points(1,:)-62.5;
Points(2,:) = Points(2,:)-62.5;

idx1 = Points(3,:)>0.01;

Points = Points(:,idx1);
Area = Area(idx1);
Normals = Normals(:,idx1);

size(Points,2)

radius = (Points(1,:).^2+Points(2,:).^2);
idx2 =(radius<(56^2)) & (Points(3,:)<59.99);

Points = Points(:,~idx2);
Area = Area(~idx2);
Normals = Normals(:,~idx2);

size(Points,2)
end


function [Points, Normals] = rotateX(Points, Normals,thetax)

Rx = [1, 0, 0;
    0,cos(thetax), -sin(thetax);
    0,sin(thetax), cos(thetax);];
Points = Rx * Points;
Normals = Rx * Normals;
end

function [Points, Normals] = rotateZ(Points, Normals,thetaz)

Rz = [cos(thetaz), -sin(thetaz), 0;
    sin(thetaz),cos(thetaz), 0;
    0,0, 1;];
Points = Rz * Points;
Normals = Rz * Normals;
end




