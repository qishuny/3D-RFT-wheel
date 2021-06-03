
%% RFT Modeling
% Qishun Yu
% 09/26/2020

% plot_alphaMap()

%% EXAMPLE START

% file_name = 'smooth_wheel_150.DXF';
% 
% file_name = 'Rover_wheel_printable_straight.DXF';
% pos = read_dxf(file_name);
% pos = pos./10;
pos = makecircle(35);
[beta,area] = initialize_beta(pos);

% plot_alphaMap()


% angular velocity
w = 0.2*pi;

T = 5;
dt = 0.1;
depth = - 2.0;
rho = 4; %????
center = [0,0];



[ForceX,ForceZ,pos_net,beta_net,area_net,alphaX_net,alphaZ_net] = netforce(depth,pos,beta,area,center);

ForceZ

figure
rectangle('Position',[-10,-6,20,4],'FaceColor','y','edgecolor','y')
scatter(pos(:,1),pos(:,2),1,'blue')
hold on
plot(center(1),center(2),'r*')
quiver(pos_net(:,1),pos_net(:,2),alphaX_net,alphaZ_net,'r')
xlim([-10,20])
ylim([-6,8])


daspect([1 1 1])
pause(0.1)
            
% simulation(T,dt,w,initial_pos,initial_beta,initial_area,initial_depth)
%EXAMPLE END
%% Functions



    function [ForceX,ForceZ,pos_net,beta_net,area_net,alphaX_net,alphaZ_net] = netforce(depth,pos_points,beta,area,center)
        pos_net = [];
        beta_net = [];
        area_net = [];
        depth_net = [];
        for i = 1:size(pos_points,1)
            if pos_points(i,2)<depth
                pos_net = [pos_net;pos_points(i,:)];
                beta_net = [beta_net;beta(i)];
                area_net = [area_net;area(i)];
                depth_net = [depth_net;abs(pos_points(i,2)-depth)];
            end
        end
        w = 0.2*pi;
        
        vel_net = vel_func(pos_net,w,center);
   
        
        figure 
        quiver(pos_net(:,1),pos_net(:,2),vel_net(:,1),vel_net(:,2))
        daspect([1 1 1])
       
%        
%         x = pos_net(:,1)';
%         y = pos_net(:,2)';
%         text1 = beta_net';
%         size(x)
%         size(y)
%         size(text1)
%         figure 
%         text(x(1,:),y(1,:),string(text1(1,:)))
%         daspect([1 1 1])
        
        
        [alphaX_net,alphaZ_net] = alpha_all(pos_net,vel_net,beta_net);
        
        ForceX = sum(alphaX_net.*area_net.*depth_net);
        ForceZ = sum(alphaZ_net.*area_net.*depth_net);
    end

% Solve for velocity at each point
    function [vel_points] = vel_func(pos_points,w,center)
        vel_points = zeros(size(pos_points));
        for iter6 = 1:size(pos_points,1)
            xtemp = pos_points(iter6,1)-center(1);
            ytemp = pos_points(iter6,2)-center(2);
            r = sqrt(xtemp^2 + ytemp^2);
            angle = atan2(ytemp,xtemp);
            vel_points(iter6,1) = r*w+(sin(angle)*r*w);
            vel_points(iter6,2) = -(cos(angle)*r*w);
        end
    end
% translation and rotation
    function [pos_new,center_new,beta] = rotate_update(pos_points,center,beta,theta,xt,zt)
        pos_new = zeros(size(pos_points));
        vel_new = zeros(size(pos_points));
        rel_points = zeros(size(pos_points));
        
        
        rel_points(:,1) = pos_points(:,1)- repmat(center(1),size(pos_points,1),1);
        rel_points(:,2) = pos_points(:,2)- repmat(center(2),size(pos_points,1),1);
        
        center_new = [center(1)+xt,center(2)+zt];
        w = 0.2*pi;
        for iter7 = 1:size(rel_points,1)
            T = [cos(theta), -sin(theta),0;
                sin(theta), cos(theta), 0;
                0,0,1];
            pTemp = T*[rel_points(iter7,1);rel_points(iter7,2);1];
            
            pos_new(iter7,:) = [pTemp(1)+center_new(1),pTemp(2)+center_new(2)];
            beta(iter7) = beta(iter7)-theta;
            if beta(iter7) >pi
                beta(iter7) = beta(iter7)-pi;
            end
        end
        
    end

% find positions of points below depth
    function [pos_out,beta_out,area_out] = pos_limit(pos_points,area,depth,beta)
        pos_out = [];
        beta_out = [];
        area_out = [];
        for i = 1:size(pos_points,1)
            if pos_points(i,2)<depth
                pos_out = [pos_out;pos_points(i,:)];
                beta_out = [beta_out;beta(i)];
                area_out = [area_out;area(i)];
            end
        end
    end

% find beta, gamma
    function [beta, gamma] = bg_func(pos_points,vel_points,beta)
        
        gamma = zeros(size(pos_points,1),1);
        for iter1 = 1:size(pos_points,1)
            
            gTemp = pi+atan2(vel_points(iter1,2),vel_points(iter1,1));
            gTemp = wrapToPi(gTemp);
            gamma(iter1) = gTemp;
        end
    end

% find total alphax and alphaz with give gamma and beta
    function [alphaX_all,alphaZ_all] = alpha_all(pos_points,vel_points,beta)
        alphaX_all = zeros(size(pos_points,1),1);
        alphaZ_all = zeros(size(pos_points,1),1);
        [beta, gamma] = bg_func(pos_points,vel_points,beta);
        
        
        beta_plot(pos_points,vel_points,beta,gamma);
        for iter2 = 1:size(pos_points,1)
            [alphaX_all(iter2),alphaZ_all(iter2)]= alpha_func(beta(iter2),gamma(iter2));
        end
    end

% find the local alphax and alphaz with give gamma and beta
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

% Read dxf file of the wheel and plot the wheel shape
    function pos_points = read_dxf(file_name)
        dxf = DXFtool(file_name);
        h = findobj(gca,'Type','line');
        x=get(h,'Xdata');
        y=get(h,'Ydata');
        
        x1 = [];
        y1 = [];
        for iter3 = 1:size(x,1)
            if iscell(x(iter3))
                xtemp = cell2mat(x(iter3));
                ytemp = cell2mat(y(iter3));
            else
                xtemp = x;
                ytemp = y;
            end
            si = size(xtemp,2);
            if si ==50
                x1 = [x1;xtemp];
                y1 = [y1;ytemp];
            elseif si ==2
                [xList, yList]= line_func (xtemp,ytemp,50);
                x1 = [x1;xList'];
                y1 = [y1;yList'];
            end
        end
        x1 = reshape(x1,[],1);
        y1 = reshape(y1,[],1);
        
        pos_points = [x1,y1];
        function [xList, yList]= line_func (x2,y2,num)
            xList = zeros(num,1);
            yList = zeros(num,1);
            for j = 1:50
                xList(j) = x2(1)+j*(x2(2)-x2(1))/(num+1);
                yList(j) = y2(1)+j*(y2(2)-y2(1))/(num+1);
            end
            
        end
    end

% DRAW A CIRCLE :)
    function pos_points = makecircle(rho)
        n = 250;
        pos_points = zeros(n,2);
        for iter4 = 1:n
            theta = 2*pi/n*iter4;
            [pos_points(iter4,1),pos_points(iter4,2)] = pol2cart(theta,rho);
        end
    end

    function [beta,area] = initialize_beta(pos_points)
        pos_copy = pos_points;
        
        [M,~] = size(pos_points);
        
        beta = zeros(size(pos_points,1),1);
        area = zeros(size(pos_points,1),1);
        for iter5  = 1:M
            temp = repmat([pos_points(iter5,1),pos_points(iter5,2)],M,1);
            a = pos_copy-temp;
            all_dist= vecnorm(a,2,2);
            [sorted_dist,sorted_index]=sort(all_dist,'ascend');
            prev_points = pos_points(sorted_index(2),:);
            next_points = pos_points(sorted_index(3),:);
            Atemp = sqrt((prev_points(1)-next_points(1))^2+(prev_points(2)-next_points(2))^2)/2;
            bTemp = pi-atan2((next_points(2)-prev_points(2)),(next_points(1)-prev_points(1)));
            if bTemp >pi
                bTemp = bTemp-pi;
            end
            beta(iter5) = bTemp;
            area(iter5) = Atemp;
        end
    end





%% Plot functions
%Plot velocity and position profile functions
    function plt_vel(pos_points,vel_points)
        figure()
        hold on
        quiver(pos_points(:,1),pos_points(:,2),vel_points(:,1),vel_points(:,2),'r')
        scatter(pos_points(:,1),pos_points(:,2),1,'blue')
        title('position and velocity')
        axis equal
        hold off
    end

% plot alpha x and alpha z as a function of beta and gamma
    function plot_alphaMap()
        
        beta = linspace(-pi/2,pi/2,500);
        gamma = linspace(-pi/2,pi/2,500);
        outputZ = zeros(size(beta,2),size(gamma,2));
        outputX = zeros(size(beta,2),size(gamma,2));
        for iter8 = 1: size(beta,2)
            for iter9 = 1:size(gamma,2)
                [alphaX, alphaZ] = alpha_func(beta(iter8),gamma(iter9));
                outputX(iter8,iter9)= alphaX;
                outputZ(iter8,iter9)= alphaZ;
            end
        end
        
        figure()
        hold on
        imagesc([beta(1,1) beta(1,end)],[gamma(1,1) gamma(1,end)],outputX);
        axis xy;
        colormap();
        colorbar;
        xlim([-pi/2,pi/2])
        ylim([-pi/2,pi/2])
        set(gca,'XTick',-0.5*pi:0.5*pi:0.5*pi)
        set(gca,'XTickLabel',{'-pi/2','0','pi/2'})
        set(gca,'YTick',-0.5*pi:0.5*pi:0.5*pi)
        set(gca,'YTickLabel',{'-pi/2','0','pi/2'})
        xlabel('\bf \gamma');
        ylabel('\bf \beta');
        title('\alpha x')
        
        hold off
        
        figure()
        hold on
        imagesc([beta(1,1) beta(1,end)],[gamma(1,1) gamma(1,end)],outputZ);
        axis xy;
        colormap();
        colorbar;
        xlim([-pi/2,pi/2])
        ylim([-pi/2,pi/2])
        xlabel('\bf\gamma');
        ylabel('\bf\beta');
        set(gca,'XTick',-0.5*pi:0.5*pi:0.5*pi)
        set(gca,'XTickLabel',{'-pi/2','0','pi/2'})
        set(gca,'YTick',-0.5*pi:0.5*pi:0.5*pi)
        set(gca,'YTickLabel',{'-pi/2','0','pi/2'})
        title('\alpha z')
        hold off
    end

% plot beta and gamma(Test)
    function beta_plot(pos_points,vel_points,beta,gamma)
        a = size(pos_points(:,1),1);
        figure()
        scatter(pos_points(:,1),pos_points(:,2),'.')
        hold on
        for iter10=1:a
            b = beta(iter10)*180/pi;
            g = gamma(iter10)*180/pi;
            text(pos_points(iter10,1),pos_points(iter10,2),num2str(b), 'FontSize', 17)
        end
%         quiver(pos_points(:,1),pos_points(:,2),vel_points(:,1),vel_points(:,2),'r')
        axis equal
    end
