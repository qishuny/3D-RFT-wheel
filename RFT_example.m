%% try to recreate 2019 Agarwal, Kamrin calculations (Fig 11); test compute_RFT_wheel output

% close all; 
% clear; clc;
r_wheel = 0.13;
w_wheel = 1.23*r_wheel;
rftwheel = [r_wheel, w_wheel];

Fg = [80, 130, 150, 190];
v0 = [1; 0; 0];

MSscale = 2.02; MMSscale = 3.05;
coeff_generic = [0.206, 0.169, 0.212, 0.358, 0.055, -0.124, 0.253, 0.007, 0.088]; % from supplementary section of original paper
Mrft = MMSscale * coeff_generic;

ang_vel = .5;
s = [-0.7, -0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5, 0.7];

ng = numel(Fg); ns = numel(s);
sinkage = zeros(ng, ns);
drawbar = zeros(ng, ns);
torques = zeros(ng, ns);

% Colors = [[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];...
%     [0.4940, 0.1840, 0.5560];[0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330];[0.6350, 0.0780, 0.1840]];

figure;
sub(1) = subplot(3,1,1); hold on;
sub(2) = subplot(3,1,2); hold on;
sub(3) = subplot(3,1,3); hold on;

for i = 1:ng
    Fgi = Fg(i);

    for j = 1:ns
        slip = s(j);
        v = (1-slip)*ang_vel*r_wheel * v0;
        
        % calculate sinkage: Fz = Fg
        fun = @(r_z) (-Fgi + dot([0,1,0], compute_RFT_wheel(rftwheel, [0;0;r_z], v, ang_vel, Mrft)));
        z_sink = fsolve(fun, -r_wheel/4);

        rft_forces = compute_RFT_wheel( rftwheel, [0;0;z_sink], v, ang_vel, Mrft );
        sinkage(i,j) = -z_sink;
        drawbar(i,j) = rft_forces(1);
        torques(i,j) = -rft_forces(3);
    end
    
    subplot(3,1,1); plot(s,drawbar(i,:));
    subplot(3,1,2); plot(s,torques(i,:)); ylabel('drive torque [Nm]')
    subplot(3,1,3); plot(s,1000*sinkage(i,:)); ylabel('sinkage [mm]')
end

linkaxes(sub,'x');
sub(3); ylim([0 inf])

subplot(3,1,1); ylabel('drawbar pull [N]'); legend(num2cell(string(Fg)),'Location','southeast');

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



