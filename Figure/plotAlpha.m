beta = linspace(-pi, pi, 1000);
gamma = linspace(-pi, pi, 1000);

[betaList,gammaList] = ndgrid(beta, gamma);
[alphaX, alphaZ] = calc_rft_alpha(betaList, gammaList, 1);


figure
tlo = tiledlayout(2,1);
h(1) = nexttile(tlo);   

mesh(gammaList,betaList,alphaX)
daspect([1 1 1])
xlabel('\bf \gamma');
ylabel('\bf \beta');
xlim([-pi pi])
xticks([-pi, -pi/2, 0,pi/2,pi])
xticklabels({'-\pi', '-\pi/2', '0','\pi/2','\pi'})
ylim([-pi pi])
yticks([-pi, -pi/2, 0,pi/2,pi])
yticklabels({'-\pi', '-\pi/2', '0','\pi/2','\pi'})

title('\bf \alpha_{x}')


h(2) = nexttile(tlo); 
mesh(gammaList,betaList,alphaZ)
daspect([1 1 1])
xlabel('\bf \gamma');
ylabel('\bf \beta');
xlim([-pi pi])
xticks([-pi, -pi/2, 0,pi/2,pi])
xticklabels({'-\pi', '-\pi/2', '0','\pi/2','\pi'})
ylim([-pi pi])
yticks([-pi, -pi/2, 0,pi/2,pi])
yticklabels({'-\pi', '-\pi/2', '0','\pi/2','\pi'})

title('\bf  \alpha_{z}')
set(h, 'Colormap', parula, 'CLim', [-1 1])
cbh = colorbar(h(end)); 
cbh.Layout.Tile = 'east'; 




function [alphaX, alphaZ] = calc_rft_alpha(beta, gamma, sf)
% using discrete Fourier transform fitting function [Li et al., 2013]
% Fourier coefficients M
M = [0.206, 0.169, 0.212, 0.358, 0.055, -0.124, 0.253, 0.007, 0.088];

alphaZ = sf .* (M(1) .* cos(0) ...
    + M(2) .* cos(2 .* beta)...
    + M(3) .* sin(2 .* beta + gamma)...
    + M(4) .* sin(gamma)...
    + M(5) .* sin((-2 .* beta) + gamma));

alphaX = sf .* (M(6) .* cos(2 .* beta + gamma)...
    + M(7) .* cos(gamma)...
    + M(8) .* cos(-2 .* beta + gamma)...
    + M(9) .* sin(2 .* beta));
end
