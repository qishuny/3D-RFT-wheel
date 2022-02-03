wheel2Ddata = matfile('smooth_wheel_125_2D.mat');
wheel3Ddata = matfile('smooth_wheel_125.mat');

Fg = [5,10, 15, 20];
s = [-0.7, -0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5, 0.7];
ng = numel(Fg); ns = numel(s);

sinkage = zeros(ng, ns);
drawbar = zeros(ng, ns);
Fz = zeros(ng, ns);
sinkage1 = zeros(ng, ns);
drawbar1 = zeros(ng, ns);
Fz1 = zeros(ng, ns);

figure;
sub(1) = subplot(3,2,1); hold on;
sub(2) = subplot(3,2,2); hold on;
sub(3) = subplot(3,2,3); hold on;
sub(4) = subplot(3,2,4); hold on;
sub(5) = subplot(3,2,5); hold on;
sub(6) = subplot(3,2,6); hold on;


rwheel = 62.5;
ang_vel = 0.5;                
wheelWidth = 60;


MMSscale = 3.05;

v0 = [1; 0; 0];
for i = 1:ng
    Fgi = Fg(i);

    for j = 1:ns
        slip = s(j);
        vcenter = (1-slip) * ang_vel * rwheel * v0;
        
        fun = @(r_z) (-Fgi + dot([0,1], RFT2Dfunc(wheel2Ddata, ang_vel, vcenter, r_z, rwheel, wheelWidth, MMSscale)));
        z_sink = fsolve(fun, rwheel/4);
        forces = RFT2Dfunc(wheel2Ddata, ang_vel, vcenter, z_sink, rwheel, wheelWidth, MMSscale);
        
        slipAngle = 0;
        wr = ang_vel;
        vc = vcenter(1);
        fun1 = @(r_z) (-Fgi + dot([0,0,1], RFT3Dfunc(wheel3Ddata, rwheel, slipAngle, wr, vc, r_z, MMSscale)));
        z_sink1 = fsolve(fun1, rwheel/4);
%         
        forces1 = RFT3Dfunc(wheel3Ddata, rwheel, slipAngle, wr, vc, z_sink1, MMSscale);
        
        sinkage1(i,j) = z_sink1;
        drawbar1(i,j) = forces1(2);
        Fz1(i,j) = forces1(3);
        sinkage(i,j) = z_sink;
        drawbar(i,j) = forces(1);
        Fz(i,j) = forces(2);

    end
    
    subplot(3,2,1); plot(s,sinkage(i,:)); ylabel('depth 2D[mm]');
    subplot(3,2,3); plot(s,drawbar(i,:)); ylabel('drawbar 2D[N]');
    subplot(3,2,5); plot(s,Fz(i,:)); ylabel('Fz 2D[N]');
    subplot(3,2,2); plot(s,sinkage1(i,:)); ylabel('depth 3D [mm]')
    subplot(3,2,4); plot(s,drawbar1(i,:)); ylabel('drawbar 3D [N]');
    subplot(3,2,6); plot(s,Fz1(i,:)); ylabel('Fz 3D [N]')

end

