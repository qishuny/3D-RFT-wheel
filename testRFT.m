clear all
close all
clc

[x1, z1] = rft_alpha(-pi/2, 0)
[x2, z2] = rft_alpha(0, 0)
[x1o, z1o] = alpha_func(-pi/2, 0)
[x2o, z2o] = alpha_func(0, 0)
[f1, f23] = normalScale(pi/4)


plot_map_original()

tic
beta = linspace(-pi, pi, 1000);
gamma = linspace(-pi, pi, 1000);

[betaList,gammaList] = ndgrid(beta, gamma);
[alphaX, alphaZ] = rft_alpha(betaList,gammaList);
toc

figure
mesh(gammaList,betaList,alphaX)
daspect([1 1 1])
xlabel('\bf \gamma');
ylabel('\bf \beta');
title('\alpha x')
colormap();
colorbar;

figure
mesh(gammaList,betaList,alphaZ)
daspect([1 1 1])
xlabel('\bf \gamma');
ylabel('\bf \beta');
title('\alpha x')
colormap();
colorbar;

function [alphaX, alphaZ] = rft_alpha(beta,gamma)
% using discrete Fourier transform fitting function [Li et al., 2013]
% beta [-pi,pi]
% gamma [-pi,pi]
% Fourier coefficients M
% granular medium: generic coefficient
% define (-1 as M1 for simplicity)
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
% scaling factor
sf = 1;
beta = wrapToPi(beta);
gamma = wrapToPi(gamma);

idxb1 = (beta >= -pi & beta <= -pi/2);
idxb2 = (beta >= pi/2 & beta <= pi);

beta(idxb1) = beta(idxb1) + pi;
beta(idxb2) = beta(idxb2) - pi;
beta = wrapToPi(beta);


idxg1 = (gamma >= -pi & gamma <= -pi/2);
idxg2 = (gamma >= pi/2 & gamma <= pi);

beta(idxg1) = -beta(idxg1);
beta(idxg2) = -beta(idxg2);

gamma(idxg1) = -pi - gamma(idxg1);
gamma(idxg2) = pi - gamma(idxg2);

gamma = wrapToPi(gamma);

alphaZ = sf .* (M(1) .* cos(0) ...
    + M(2) .* cos(2 .* beta)...
    + M(3) .* sin(2 .* beta + gamma)...
    + M(4) .* sin(gamma)...
    + M(5) .* sin((-2 .* beta) + gamma));
alphaX = sf .* (M(6) .* cos(2 .* beta + gamma)...
    + M(7) .* cos(gamma)...
    + M(8) .* sin(-2 .* beta + gamma)...
    + M(9) .* sin(2 .* beta));

alphaX(idxg1) = -alphaX(idxg1);
alphaX(idxg2) = -alphaX(idxg2);
end


function [f1, f23] = normalScale(phi)
a1 = 0.44;
a2 = 3.62;
a3 = 1.61;
a4 = 0.41;

b1 = 1.99;
b2 = 1.61;
b3 = 0.97;
b4 = 4.31;

phi = wrapToPi(phi);

if phi>pi/2 && phi<=pi
    phi = pi-phi;
elseif phi<0 && phi>=-pi/2
    phi = -phi;
elseif phi<-pi/2 && phi >=-pi
    phi = phi+pi;
end
v1v = sin(phi);
v23v = cos(phi);

F1 = a1 * tanh(a2 * v1v - a3) + a4;
F23 = b1 * atanh(b2 * v23v - b3) + b4;
F1max =a1 * tanh(a2 * sin(pi/2) - a3) + a4;
F23max = b1 * atanh(b2 * cos(0) - b3) + b4;

f1 = F1 / F1max;
f23 = F23 / F23max;
end

% find the local alphax and alphaz with give gamma and beta
function [alphaX, alphaZ] = alpha_func(beta,gamma)
% using discrete Fourier transform fitting function
% beta [-pi,pi]
% gamma [-pi,pi]
% Fourier coefficients M
% granular medium: generic coefficient
% define (-1 as M1 for simplicity)
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
% scaling factor
sf = 1;


if beta >= -pi && beta <= -pi/2
    beta = beta + pi;
elseif beta >= pi/2 && beta <= pi
    beta = beta - pi;
end

if gamma >= -pi && gamma <= -pi/2
    alphaZ = sf*(M(1)*cos(0)...
        +M(2)*cos(2*(-beta))...
        +M(3)*sin(2*(-beta)+(-pi-gamma))...
        +M(4)*sin((-pi-gamma))+M(5)*sin((-2*(-beta))+(-pi-gamma)));
    alphaX = -sf*(M(6)*cos(2*(-beta)+(-pi-gamma))+M(7)*cos((-pi-gamma))...
        +M(8)*sin(-2*(-beta)+(-pi-gamma))+M(9)*sin(2*(-beta)));
elseif gamma >= pi/2 && gamma <= pi
    alphaZ = sf*(M(1)*cos(0)...
        +M(2)*cos(2*(-beta))...
        +M(3)*sin(2*(-beta)+(pi-gamma))...
        +M(4)*sin((pi-gamma))...
        +M(5)*sin((-2*(-beta))+(pi-gamma)));
    alphaX = -sf*(M(6)*cos(2*(-beta)+(pi-gamma))+M(7)*cos((pi-gamma))...
        +M(8)*sin(-2*(-beta)+(pi-gamma))+M(9)*sin(2*(-beta)));
else
    alphaZ = sf*(M(1)*cos(0)...
        +M(2)*cos(2*beta)...
        +M(3)*sin(2*beta+gamma)...
        +M(4)*sin(gamma)...
        +M(5)*sin((-2*beta)+gamma));
    alphaX = sf*(M(6)*cos(2*beta+gamma)+M(7)*cos(gamma)...
        +M(8)*sin(-2*beta+gamma)+M(9)*sin(2*beta));
end

end

function plot_map_original()

beta = linspace(-pi,pi,1000);
gamma = linspace(-pi,pi,1000);
outputZ = zeros(size(beta,2),size(gamma,2));
outputX = zeros(size(beta,2),size(gamma,2));
tic
for iter8 = 1: size(beta,2)
    for iter9 = 1:size(gamma,2)
        [alphaX, alphaZ] = alpha_func(beta(iter8),gamma(iter9));
        outputX(iter8,iter9)= alphaX;
        outputZ(iter8,iter9)= alphaZ;
    end
end
toc
figure()
hold on
imagesc([beta(1,1) beta(1,end)],[gamma(1,1) gamma(1,end)],outputX);
axis xy;
colormap();
colorbar;
xlim([-pi,pi])
ylim([-pi,pi])
set(gca,'XTick',-0.5*pi:0.5*pi:0.5*pi)
set(gca,'XTickLabel',{'-pi/2','0','pi/2'})
set(gca,'YTick',-0.5*pi:0.5*pi:0.5*pi)
set(gca,'YTickLabel',{'-pi/2','0','pi/2'})
xlabel('\bf \gamma');
ylabel('\bf \beta');
title('\alpha x')

hold off

figure()
hold on
imagesc([beta(1,1) beta(1,end)],[gamma(1,1) gamma(1,end)],outputZ);
axis xy;
colormap();
colorbar;
xlim([-pi,pi])
ylim([-pi,pi])
xlabel('\bf\gamma');
ylabel('\bf\beta');
set(gca,'XTick',-0.5*pi:0.5*pi:0.5*pi)
set(gca,'XTickLabel',{'-pi/2','0','pi/2'})
set(gca,'YTick',-0.5*pi:0.5*pi:0.5*pi)
set(gca,'YTickLabel',{'-pi/2','0','pi/2'})
title('\alpha z')
hold off
end