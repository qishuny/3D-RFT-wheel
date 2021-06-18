[x1, z1] = rft_alpha(pi/2, 0)
[x2, z2] = rft_alpha(pi/2, pi)

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
if beta >= -pi && beta <= -pi/2
    beta = beta + pi;
elseif beta >= pi/2 && beta <= pi
    beta = beta - pi;
end
beta = wrapToPi(beta);
gamma = wrapToPi(gamma);
if gamma >= -pi && gamma <= -pi/2
    alphaZ = sf*(M(1)*cos(0)+M(2)*cos(2*(-beta))+M(3)*sin(2*(-beta)+(-pi-gamma))...
        +M(4)*sin((-pi-gamma))+M(5)*sin((-2*(-beta))+(-pi-gamma)));
    alphaX = -sf*(M(6)*cos(2*(-beta)+(-pi-gamma))+M(7)*cos((-pi-gamma))...
        +M(8)*sin(-2*(-beta)+(-pi-gamma))+M(9)*sin(2*(-beta)));
elseif gamma >= pi/2 && gamma <= pi
    alphaZ = sf*(M(1)*cos(0)+M(2)*cos(2*(-beta))+M(3)*sin(2*(-beta)+(pi-gamma))...
        +M(4)*sin((pi-gamma))+M(5)*sin((-2*(-beta))+(pi-gamma)));
    alphaX = -sf*(M(6)*cos(2*(-beta)+(pi-gamma))+M(7)*cos((pi-gamma))...
        +M(8)*sin(-2*(-beta)+(pi-gamma))+M(9)*sin(2*(-beta)));
else
    alphaZ = sf*(M(1)*cos(0)+M(2)*cos(2*beta)+M(3)*sin(2*beta+gamma)...
        +M(4)*sin(gamma)+M(5)*sin((-2*beta)+gamma));
    alphaX = sf*(M(6)*cos(2*beta+gamma)+M(7)*cos(gamma)...
        +M(8)*sin(-2*beta+gamma)+M(9)*sin(2*beta));
end

end