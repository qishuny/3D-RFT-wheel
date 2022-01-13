clear all
close all
clc

wheeldata = matfile('smooth_wheel_130_2D.mat');
s = [-0.7, -0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5, 0.7];
wr = 0.5;

ang_vel = 0.5;
r_wheel = 65;
sinkage = 27.4;
slip = s(1);
vcenter = (1-slip) * ang_vel * r_wheel * 1;
wheelWidth = 1.23 * 130;
MMSscale = 3.05;

[Fx, Fz] = RFT2Dfunc(wheeldata, wr, vcenter, sinkage, wheelWidth, MMSscale)
