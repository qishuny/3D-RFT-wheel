wheel3Ddata = matfile('smooth_wheel_125.mat');


slipAngle = 45 * pi / 180;
wr = 0;
vc = 10;
MMSscale = 0.6;
z_sink1 = 50;

forces1 = RFT3Dfunc(wheel3Ddata, 62.5, slipAngle, wr, vc, z_sink1, MMSscale)
        