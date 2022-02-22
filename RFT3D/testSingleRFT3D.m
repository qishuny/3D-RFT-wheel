wheel3Ddatas = matfile('smooth_wheel_125.mat');
wheel3Ddatag = matfile('grousered_wheel_125.mat');
slipAngle = 45 * pi / 180;
wr = 0;
vc = 10;
MMSscale = 0.6;
z_sink1 = 20;

forces1 = RFT3DNewfunc(wheel3Ddatag, 62.5, slipAngle, wr, vc, z_sink1, MMSscale, 0);
% forces = RFT3Dfunc(wheel3Ddata, 62.5, slipAngle, wr, vc, z_sink1, MMSscale) 
[forces] = RFT3DDEMfunc(wheel3Ddatag, 62.5, slipAngle, wr, vc, z_sink1, MMSscale, 1);

% pointListg = wheel3Ddata.Points
% [idxOut, depthList, pile, under] = SandDeformation(pointList, slipAngle, depth)