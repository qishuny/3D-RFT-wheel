clear all
close all
clc

% read element information generated from meshstuff.m
wheel = matfile('wheel.mat');
area = matfile('area.mat');
normal = matfile('normal.mat');

% element points information (3*n) [x;y;z]
pointList = wheel.pointList;
% element areas information (1*n)
areaList = area.areaList;
% element normal vectors information (3*n)
normalList = normal.normalList;

% element radius to the center information (1*n)
rList = sqrt(pointList(2,:).^2+pointList(3,:).^2);
% element tangent angle (1*n)
angleList = atan2(pointList(3,:),pointList(2,:))+pi/2;


% angular speed mm/s
w = 2;

% the velocity of the center of the wheel mm/s
vcorx = 200;
vcory = 0;
vcorz = 0;
vx = zeros([1,size(angleList,2)])+vcorx;
vz = sin(angleList).*rList.*w+vcory;
vy = cos(angleList).*rList.*w+vcorz;



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