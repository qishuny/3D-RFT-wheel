depth = [50.0267   32.0452   23.3210   14.7238   14.0679   17.5152   22.3674   25.0873   25.3562];
% depth = 4 * [10 10 10 10 10 10 10 10 10];
Vry = [0.00001 3 5 8 10 12 20 40 100];
slip = [-1 -0.7 -0.5 -0.2 0 0.23 0.5 0.75 0.9];
scale = 0.6985;
wheeldata = matfile('../data/smooth_wheel_125.mat');
wheel2Ddata = matfile('smooth_wheel_125_2D.mat');
r_wheel = 0.0625;
w_wheel = 0.06;
rftwheel = [r_wheel, w_wheel];
vs = [0.01; 0; 0];
coeff_generic = [0.206, 0.169, 0.212, 0.358, 0.055, -0.124, 0.253, 0.007, 0.088]; % from supplementary section of original paper
Mrft = scale * coeff_generic;

drawbar = zeros(9, 1);
Fz = zeros(9, 1);
drawbar1 = zeros(9, 1);
Fz1 = zeros(9, 1);
drawbar2 = zeros(9, 1);
Fz2 = zeros(9, 1);
for i = 1:9
    ang_vel = Vry(i)/62.5;
    z_sink1 = -depth(i)/1000;
    forces = compute_RFT_wheel(rftwheel, [0;0;z_sink1], vs, ang_vel, Mrft )
    drawbar(i) = forces(1);
    Fz(i) = forces(2);
    
    radius = 62.5;
    vcenter = 10;
    slipAngle = 0;
    sinkage = depth(i);
    forcenew = RFT3DNewfunc(wheeldata, radius, slipAngle, ang_vel, vcenter, sinkage, scale);
    [Force] = RFT3Dfunc(wheeldata, radius, slipAngle, ang_vel, vcenter, sinkage, scale);
    drawbar1(i) = Force(2);
    Fz1(i) = Force(3);
    drawbar2(i) = forcenew(2);
    Fz2(i) = forcenew(3);
    
end
figure;

subplot(3,2,1); plot(slip,drawbar); ylabel('Fx[mm]');
subplot(3,2,2); plot(slip,Fz); ylabel('Fz[mm]');
subplot(3,2,3); plot(slip,drawbar1); ylabel('Fx 3D[mm]');
subplot(3,2,4); plot(slip,Fz1); ylabel('Fz 3D[mm]');
subplot(3,2,5); plot(slip,drawbar2); ylabel('Fx new 3D[mm]');
subplot(3,2,6); plot(slip,Fz2); ylabel('Fz new 3D[mm]');

figure;
plot(slip, drawbar);
hold on
plot(slip, drawbar2);
legend("2D", "3D");
figure;
plot(slip, Fz);
hold on
plot(slip, Fz2);
legend("2D", "3D");
%%
function forces = compute_RFT_wheel( wheel_geometry, pos, vel, angvel, M_rft )
%COMPUTE_FT_WHEEL compute drag forces & torques associated with dynamic wheel
% wheel is located in vehicle coordinates (+x to front, +z up)
% wheel is free to move in 2D, and driving torque is applied in y direction
% pos & vel are 3D vectors; ang_vel is a single value about +y
 
    r_wheel = wheel_geometry(1); w_wheel = wheel_geometry(2);
    sink = -pos(3);
    
    % z=0 where the wheel axle is exactly the wheel's radius off the ground, 
    % i.e. the wheel is resting exactly on the ground surface with no sinkage
    if sink < 1e-5
%         disp('There is no ground contact')
        forces = zeros(1,3);
        return
    end
    
        
    ang_step = pi/30; % dA = (r_wheel * ang_step) * w_wheel;
    if sink >= 2*r_wheel % wheel is fully sunk; RFT does not necessarily operate here
%         theta_lead = 3*pi/2;
        theta_lead = NaN;
    else % limit sampling to area in contact with sand
        theta_lead = asin(1 - sink/r_wheel);
    end
    if norm(vel) > 1e-3
        theta_vel_bound = atan2(vel(1), vel(3));
        if theta_vel_bound < 0
            theta_vel_bound = theta_vel_bound + 2*pi;
        end
    else
        theta_vel_bound = NaN;
    end
    
    if isnan(theta_lead) && isnan(theta_vel_bound)
        theta_min = 0; theta_max = 2*pi;
    else
        theta_max = min([theta_vel_bound, pi-theta_lead]);
        theta_min = max([theta_vel_bound-pi, theta_lead]);
        if theta_min >= theta_max
%             disp('no leading edge contact')
            forces = zeros(1,3);
            return
        end
    end

    % centers test points about theta = pi/2 (straight down) and with fixed step size
    theta_range_full = [flip(pi/2:-ang_step:0), (pi/2+ang_step):ang_step:2*pi-ang_step];
    theta_range = theta_range_full(theta_range_full > theta_min & theta_range_full < theta_max);
    theta_range = [theta_min, theta_range, theta_max];
    
    it = size(theta_range,2);
      
    
    % set up RFT calc
    Torq = [0;0;0];
    Fx = 0;
    Fz = 0;
    omega = [0; angvel; 0];
    
    for i = 1:it-1
        theta = (theta_range(i) + theta_range(i+1))/2;
        dA = (r_wheel * (theta_range(i+1)-theta_range(i))) * w_wheel;
        
        ri = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)] * [r_wheel;0;0];
        vi = vel + cross(omega, ri);
        
        z_sink = ri(3) + pos(3) + r_wheel;
        
        if norm(vi) >= 5e-4 && dot(ri, vi) >= -1e-5 && z_sink < 0
            gamma = atan2(-vi(3), -vi(1));
            beta = theta + pi/2;
%             disp([beta, gamma])

            [dFx, dFz] = computeRF( beta, gamma, M_rft, dA*10^4, -z_sink*100 );
%             [ax, az] = computeAlphasMat( beta, gamma, M_rft );
%             disp([ax, az])
%             disp([dA, z_sink])
            
            Fx = Fx + dFx; Fz = Fz + dFz;
            Torq = Torq + cross(ri, [dFx; 0; +dFz]); 
%             disp([ri, [dFx;0;dFz]])
        end

    end
    Torq = Torq(2);
    forces = [Fx, Fz, Torq];  
    
  
end

function [dFx, dFz] = computeRF( beta, gamma, M, dA,z)
    %% Calculate alpha z 
    alphaZ = M(5).*sin(2*pi*((-beta./pi)+(gamma/(2*pi))));      % m = -1, n = 1
    alphaZ = alphaZ + M(1).*cos(2*pi*0);                        % m = 0,  n = 0
    alphaZ = alphaZ + M(4).*sin(2*pi*(gamma/(2*pi)));             % m = 0,  n = 1
    alphaZ = alphaZ + M(2).*cos(2*pi*(beta./pi));               % m = 1,  n = 0
    alphaZ = alphaZ + M(3).*sin(2*pi*((beta/pi)+(gamma/(2*pi)))); % m = 1,  n = 1
    
    %% Calculate alpha x
    alphaX = M(8).*cos(2*pi*((-beta./pi)+(gamma/(2.*pi))));     % m = -1, n = 1    
    alphaX = alphaX + M(7).*cos(2*pi*(gamma/(2.*pi)));          % m = 0, n = 1 
    alphaX = alphaX + M(9).*sin(2*pi*(beta./pi));               % m = 1, n = 0 
    alphaX = alphaX + M(6).*cos(2*pi*((beta/pi)+(gamma/(2*pi)))); % m = 1, n = 1
    
    %% Calculate forces
    dFz = z*dA*alphaZ;
    dFx = z*dA*alphaX;

end
