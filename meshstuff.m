clear all
close all
clc


model = createpde();
gc = importGeometry(model,'Rover_wheel_printable_straight.stl');

figure
pdegplot(model,'FaceLabels','on')
rotate(gc,180,[0 53 0],[53 0 0]);
rotate(gc,180,[0 0 0],[50 0 0]);
rotate(gc,90,[0 0 0],[0 1 0]);

mesh = generateMesh(model,'Hmax',2,'GeometricOrder','linear');

figure
pdeplot3D(model)

pointList = [];
areaList = [];
normalList = [];

for i = 117:144
    [pointarray, areaarray, normalarray] = generatePointsfromFace(mesh,i);
    pointList = [pointList pointarray];
    areaList = [areaList areaarray];
    normalList = [normalList normalarray];
end
for i = 42:113
    [pointarray, areaarray, normalarray] = generatePointsfromFace(mesh,i);
    pointList = [pointList pointarray];
    areaList = [areaList areaarray];
    normalList = [normalList normalarray];
end
for i = 2:26
    [pointarray, areaarray, normalarray] = generatePointsfromFace(mesh,i);
    pointList = [pointList pointarray];
    areaList = [areaList areaarray];
    normalList = [normalList normalarray];
end

figure
plot3(pointList(1,:),pointList(2,:),pointList(3,:),'ok','MarkerFaceColor','g')
daspect([1 1 1])

function [pointarray, areaarray, normalarray] = generatePointsfromFace(mesh,faceID)
    pointarray = [];
    areaarray = [];
    normalarray = [];
    
    Nf2 = findNodes(mesh,'region','Face',faceID);
    list1 = findElements(mesh,'attached',Nf2);
    for i = 1:size(list1,2)
        Lia = ismember(mesh.Elements(:,list1(i))',Nf2);
        if sum(Lia)==3
            temp = [];
            for j = 1:4
                if(Lia(j)==1)    
                    idx = mesh.Elements(j,list1(i));       
                    temp = [temp [mesh.Nodes(1,idx);mesh.Nodes(2,idx);mesh.Nodes(3,idx)]];
                end
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
            
            %normal vector
            n = [((y2-y1)*(z3-z1)-(y3-y1)*(z2-z1));
                ((z2-z1)*(x3-x1)-(x2-x1)*(z3-z1));
                ((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));];
            
            pointarray = [pointarray [mean(temp(1,:)); mean(temp(2,:));mean(temp(3,:))]];
            areaarray = [areaarray AreaTemp];
            normalarray= [normalarray n];
        end

    end
end
