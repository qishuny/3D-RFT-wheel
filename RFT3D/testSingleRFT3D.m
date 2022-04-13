wheel3Ddatas = matfile('smooth_wheel_125.mat');
wheel3Ddatag = matfile('grousered_wheel_125.mat');
slipAngle = 22.5 * pi / 180;
wr = 0;
vc = 10;
MMSscale = 0.6;
<<<<<<< HEAD
z_sink1 = 50;
=======
z_sink1 = 20;
>>>>>>> parent of d84a5e0... retry Brian's model

% forces1 = RFT3Dfunc(wheel3Ddatag, 62.5, slipAngle, wr, vc, z_sink1, MMSscale, 0);

[forces] = RFT3DSandfunc(wheel3Ddatas, 62.5, slipAngle, wr, vc, z_sink1, MMSscale, 1);
