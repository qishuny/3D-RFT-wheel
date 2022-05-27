wheel3Ddatas = matfile('smooth_wheel_125.mat');
wheel3Ddatag = matfile('grousered_wheel_125.mat');
slipAngle = 22.5 * pi / 180;
wr = 0;
vc = 10;
MMSscale = 0.6;
z_sink1 = 15;
radius = 62.5;
w = wr / radius;
% forces1 = RFT3Dfunc(wheel3Ddatag, 62.5, slipAngle, w, vc, z_sink1, MMSscale, 0);

[forces] = RFT3DSandfunc(wheel3Ddatag, 62.5, slipAngle, w, vc, z_sink1, wr, MMSscale, 1);
