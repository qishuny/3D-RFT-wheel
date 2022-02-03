clear all
close all

% radius in mm
rho = 62.5;
[pos_points, area_points, normal_points] = makecircle(rho);
Points = pos_points';
Area = area_points';
Normals = normal_points';
save('smooth_wheel_125_2D.mat', 'Points', 'Area', 'Normals');

figure
plot(Points(1,:),Points(2,:), 'ok','MarkerFaceColor',[0,0.5,0.5])
hold on
quiver(Points(1,:),Points(2,:),Normals(1,:),Normals(2,:),1,'Color', [0,0.2,0.8]);
daspect([1 1 1])


function [pos_points, area_points, normal_points] = makecircle(rho)
n = 250;
pos_points = zeros(n,2);
area_points= (pi * 2 * rho / n) * ones(n,1);
normal_points = zeros(n,2);

for iter4 = 1:n
    theta = 2*pi/n*iter4;
    [pos_points(iter4,1),pos_points(iter4,2)] = pol2cart(theta,rho);
    normal_points(iter4,1) = cos(theta);
    normal_points(iter4,2) = sin(theta);
    
end
end


