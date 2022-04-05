wheel3Ddatas = matfile('smooth_wheel_125.mat');
wheel3Ddatag = matfile('grousered_wheel_125.mat');
slipAngle = 22 * pi / 180;
wr = 0;
vc = 10;
MMSscale = 0.6;
z_sink1 = 15;

% forces1 = RFT3Dfunc(wheel3Ddatag, 62.5, slipAngle, wr, vc, z_sink1, MMSscale, 0);

[forces] = RFT3DSandfunc(wheel3Ddatas, 62.5, slipAngle, wr, vc, z_sink1, MMSscale, 1);
