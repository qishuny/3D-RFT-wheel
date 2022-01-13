sf1 = 0.175;
sf2 = 0.5;



syms theta r_wheel ri vi vx omega angvel

vel = [vx; 0; 0];
omega = [0; angvel; 0];
ri = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)] * [r_wheel;0;0]
vi = vel + cross(omega, ri)

% beta = -pi;
% gamma = -pi;
% sf = 1;
M = [0.206, 0.169, 0.212, 0.358, 0.055, -0.124, 0.253, 0.007, 0.088];
% dA = 1;
% z = 1;
% 
% [alphaX, alphaZ] = calc_rft_alpha(beta, gamma, sf)
% [dFx, dFz] = computeRF(beta, gamma, M, dA,z)
% 
% plot_alphaMap()


wheeldata = matfile('smooth_wheel_130_2D.mat');
Fg = [80, 130, 150, 190];
s = [-0.7, -0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5, 0.7];

ang_vel = 0.5;
v0 = [1; 0; 0];

MMSscale = 3.05;
Fgi = Fg(1);
slip = s(1);
rwheel = 130;
vcenter = (1-slip) * ang_vel * rwheel * v0;
wheelWidth = 1.23 * rwheel;
% calculate sinkage: Fz = Fg
fun = @(r_z) (-Fgi + dot([0,1], RFT2Dfunc(wheeldata, ang_vel, vcenter, r_z, rwheel, wheelWidth, MMSscale)));
z_sink = fsolve(fun, rwheel/4);

forces = RFT2Dfunc(wheeldata, ang_vel, vcenter, z_sink, rwheel, wheelWidth, MMSscale);
forces(1)
forces(2)
z_sink



function [alphaX, alphaZ] = calc_rft_alpha(beta, gamma, sf)
% using discrete Fourier transform fitting function [Li et al., 2013]
% Fourier coefficients M
A00 = 0.206;
A10 = 0.169;
B11 = 0.212;
B01 = 0.358;
BM11 = 0.055;
C11 = -0.124;
C01 = 0.253;
CM11 = 0.007;
D10 = 0.088;
M = [A00, A10, B11, B01, BM11, C11, C01, CM11, D10];

% beta = wrapToPi(beta);
% gamma = wrapToPi(gamma);

% idxb1 = (beta >= -pi & beta <= -pi/2);
% idxb2 = (beta >= pi/2 & beta <= pi);
% 
% beta(idxb1) = beta(idxb1) + pi;
% beta(idxb2) = beta(idxb2) - pi;
% beta = wrapToPi(beta);
% 
% idxg1 = (gamma >= -pi & gamma <= -pi/2);
% idxg2 = (gamma >= pi/2 & gamma <= pi);
% 
% beta(idxg1) = -beta(idxg1);
% beta(idxg2) = -beta(idxg2);
% 
% gamma(idxg1) = -pi - gamma(idxg1);
% gamma(idxg2) = pi - gamma(idxg2);
% 
% gamma = wrapToPi(gamma);

alphaZ = sf .* (M(1) .* cos(0) ...
    + M(2) .* cos(2 .* beta)...
    + M(3) .* sin(2 .* beta + gamma)...
    + M(4) .* sin(gamma)...
    + M(5) .* sin((-2 .* beta) + gamma));

alphaX = sf .* (M(6) .* cos(2 .* beta + gamma)...
    + M(7) .* cos(gamma)...
    + M(8) .* cos(-2 .* beta + gamma)...
    + M(9) .* sin(2 .* beta));

% alphaX(idxg1) = -alphaX(idxg1);
% alphaX(idxg2) = -alphaX(idxg2);
end


function [dFx, dFz] = computeRF( beta, gamma, M, dA,z)
    %% Calculate alpha z 
    alphaZ = M(5).*sin(2*pi*((-beta./pi)+(gamma/(2*pi))));      % m = -1, n = 1
    alphaZ = alphaZ + M(1).*cos(2*pi*0);                        % m = 0,  n = 0
    alphaZ = alphaZ + M(4).*sin(2*pi*(gamma/(2*pi)));             % m = 0,  n = 1
    alphaZ = alphaZ + M(2).*cos(2*pi*(beta./pi));               % m = 1,  n = 0
    alphaZ = alphaZ + M(3).*sin(2*pi*((beta/pi)+(gamma/(2*pi)))); % m = 1,  n = 1
    
    %% Calculate alpha x
    alphaX = M(8).*cos(2*pi*((-beta./pi)+(gamma/(2.*pi))));     % m = -1, n = 1    
    alphaX = alphaX + M(7).*cos(2*pi*(gamma/(2.*pi)));          % m = 0, n = 1 
    alphaX = alphaX + M(9).*sin(2*pi*(beta./pi));               % m = 1, n = 0 
    alphaX = alphaX + M(6).*cos(2*pi*((beta/pi)+(gamma/(2*pi)))); % m = 1, n = 1
    
    %% Calculate forces
    dFz = z*dA*alphaZ;
    dFx = z*dA*alphaX;

end


function plot_alphaMap()

[X1,X2] = meshgrid(-pi:.2:pi);

sizeX = size(X1,1);
sizeY = size(X2,2);

alphaXS = zeros(size(X1));
alphaZS = zeros(size(X1));


sf = 1;
M = [0.206, 0.169, 0.212, 0.358, 0.055, -0.124, 0.253, 0.007, 0.088];
dA = 1;
z = 1;

for i = 1:sizeX
    for j = 1:sizeY
       [alphaXS(i,j), alphaZS(i,j)] = computeRF(X2(i,j), X1(i,j),M,dA,z); 
    end
end

[alphaX, alphaZ] = calc_rft_alpha(X2,X1,1);
figure()
mesh(X1,X2,alphaX)
xlim([-pi,pi])
ylim([-pi,pi])
colorbar;
set(gca,'XTick',-0.5*pi:0.5*pi:0.5*pi)
set(gca,'XTickLabel',{'-pi/2','0','pi/2'})
set(gca,'YTick',-0.5*pi:0.5*pi:0.5*pi)
set(gca,'YTickLabel',{'-pi/2','0','pi/2'})
xlabel('\bf \gamma');
ylabel('\bf \beta');
title('\alpha x')

figure()
mesh(X1,X2,alphaXS)
xlim([-pi,pi])
ylim([-pi,pi])
colorbar;
set(gca,'XTick',-0.5*pi:0.5*pi:0.5*pi)
set(gca,'XTickLabel',{'-pi/2','0','pi/2'})
set(gca,'YTick',-0.5*pi:0.5*pi:0.5*pi)
set(gca,'YTickLabel',{'-pi/2','0','pi/2'})
xlabel('\bf \gamma');
ylabel('\bf \beta');
title('\alpha xs')

figure()
mesh(X1,X2,alphaZ)
xlim([-pi,pi])
ylim([-pi,pi])
colorbar;
set(gca,'XTick',-0.5*pi:0.5*pi:0.5*pi)
set(gca,'XTickLabel',{'-pi/2','0','pi/2'})
set(gca,'YTick',-0.5*pi:0.5*pi:0.5*pi)
set(gca,'YTickLabel',{'-pi/2','0','pi/2'})
xlabel('\bf \gamma');
ylabel('\bf \beta');
title('\alpha y')

figure()
mesh(X1,X2,alphaZS)
xlim([-pi,pi])
ylim([-pi,pi])
colorbar;
set(gca,'XTick',-0.5*pi:0.5*pi:0.5*pi)
set(gca,'XTickLabel',{'-pi/2','0','pi/2'})
set(gca,'YTick',-0.5*pi:0.5*pi:0.5*pi)
set(gca,'YTickLabel',{'-pi/2','0','pi/2'})
xlabel('\bf \gamma');
ylabel('\bf \beta');
title('\alpha zs')

% figure()
% hold on
% imagesc([beta(1,1) beta(1,end)],[gamma(1,1) gamma(1,end)],outputX);
% axis xy;
% colormap();
% colorbar;
% xlim([-pi/2,pi/2])
% ylim([-pi/2,pi/2])
% set(gca,'XTick',-0.5*pi:0.5*pi:0.5*pi)
% set(gca,'XTickLabel',{'-pi/2','0','pi/2'})
% set(gca,'YTick',-0.5*pi:0.5*pi:0.5*pi)
% set(gca,'YTickLabel',{'-pi/2','0','pi/2'})
% xlabel('\bf \gamma');
% ylabel('\bf \beta');
% title('\alpha x')
% 
% hold off
% 
% figure()
% hold on
% imagesc([beta(1,1) beta(1,end)],[gamma(1,1) gamma(1,end)],outputZ);
% axis xy;
% colormap();
% colorbar;
% xlim([-pi/2,pi/2])
% ylim([-pi/2,pi/2])
% xlabel('\bf\gamma');
% ylabel('\bf\beta');
% set(gca,'XTick',-0.5*pi:0.5*pi:0.5*pi)
% set(gca,'XTickLabel',{'-pi/2','0','pi/2'})
% set(gca,'YTick',-0.5*pi:0.5*pi:0.5*pi)
% set(gca,'YTickLabel',{'-pi/2','0','pi/2'})
% title('\alpha z')
% hold off
end