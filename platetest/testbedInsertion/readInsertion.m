load('Zs.mat');
load('Pzs.mat');
load('expfit.mat');
% N/cm^3

coeff1 = expfit.a;
coeff2 = expfit.b;
alphaZ = Pzs ./ Zs ./ 1000;
alphaZfit = coeff1 .* (Zs .^ coeff2) ./ Zs ./ 1000;
figure
plot(Zs, Pzs);
title('Pzs as a function of Zs');

figure
plot(Zs, alphaZ);
hold on
plot(Zs, alphaZfit);
title('alphaZ as a function of Zs');