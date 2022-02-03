load('output/all_smooth_data_2.mat')
% load('output/RFToutput.mat')
% load('output/RFTDEMsameDepthoutput.mat')
% load('output/RFTDEMoutput.mat')
% load('output/RFTDEMPoly55output.mat')
% load('output/RFTbalanceoutput.mat')

load('output/RFTtestoutput.mat')
% load('output/RFTDEMoptoutput.mat')


x = [-1 -0.7 -0.5 -0.2 0 0.23 0.5 0.75 0.9];

Fx = cell(9,7);
Fy = cell(9,7);
Fz = cell(9,7);

FxRFT = cell(9,7);
FyRFT = cell(9,7);
FzRFT = cell(9,7);

depth = cell(9,7);
depthRFT = cell(9,7);

for i=1:length(all_results)
    exp_result = all_results(i);
    rft_result = RFToutput(i);
    %Select color of plotted point
    [column] = switchstuff(exp_result.Vry);
    [column_rft] = switchstuff(rft_result.wr);
    switch exp_result.beta
        case 0
                    
            Fx{column, 1} = [Fx{column, 1} -exp_result.avg_Fy];
            Fy{column, 1} = [Fy{column, 1} -exp_result.avg_Fx];
            Fz{column, 1} = [Fz{column, 1} -exp_result.avg_Fz];
            
            FxRFT{column, 1} = [FxRFT{column, 1} rft_result.ForceX];
            FyRFT{column, 1} = [FyRFT{column, 1} rft_result.ForceY];
            FzRFT{column, 1} = [FzRFT{column, 1} rft_result.ForceZ];
            
            depth{column, 1} = [depth{column, 1} -exp_result.avg_Z];
            depthRFT{column, 1} = [depthRFT{column, 1} -rft_result.depth];
        case 15
            
            Fx{column, 2} = [Fx{column, 2} -exp_result.avg_Fy];
            Fy{column, 2} = [Fy{column, 2} -exp_result.avg_Fx];
            Fz{column, 2} = [Fz{column, 2} -exp_result.avg_Fz];
            
            FxRFT{column, 2} = [FxRFT{column, 2} rft_result.ForceX];
            FyRFT{column, 2} = [FyRFT{column, 2} rft_result.ForceY];
            FzRFT{column, 2} = [FzRFT{column, 2} rft_result.ForceZ];
            
            depth{column, 2} = [depth{column, 2} -exp_result.avg_Z];
            depthRFT{column, 2} = [depthRFT{column, 2} -rft_result.depth];
            
        case 30
            
            Fx{column, 3} = [Fx{column, 3} -exp_result.avg_Fy];
            Fy{column, 3} = [Fy{column, 3} -exp_result.avg_Fx];
            Fz{column, 3} = [Fz{column, 3} -exp_result.avg_Fz];
            
            FxRFT{column, 3} = [FxRFT{column, 3} rft_result.ForceX];
            FyRFT{column, 3} = [FyRFT{column, 3} rft_result.ForceY];
            FzRFT{column, 3} = [FzRFT{column, 3} rft_result.ForceZ];
            
            depth{column, 3} = [depth{column, 3} -exp_result.avg_Z];
            depthRFT{column, 3} = [depthRFT{column, 3} -rft_result.depth];
            
        case 45
            
            Fx{column, 4} = [Fx{column, 4} -exp_result.avg_Fy];
            Fy{column, 4} = [Fy{column, 4} -exp_result.avg_Fx];
            Fz{column, 4} = [Fz{column, 4} -exp_result.avg_Fz];
            
            FxRFT{column, 4} = [FxRFT{column, 4} rft_result.ForceX];
            FyRFT{column, 4} = [FyRFT{column, 4} rft_result.ForceY];
            FzRFT{column, 4} = [FzRFT{column, 4} rft_result.ForceZ];
            
            depth{column, 4} = [depth{column, 4} -exp_result.avg_Z];
            depthRFT{column, 4} = [depthRFT{column, 4} -rft_result.depth];
            
        case 60
            
            Fx{column, 5} = [Fx{column, 5} -exp_result.avg_Fy];
            Fy{column, 5} = [Fy{column, 5} -exp_result.avg_Fx];
            Fz{column, 5} = [Fz{column, 5} -exp_result.avg_Fz];
            
            FxRFT{column, 5} = [FxRFT{column, 5} rft_result.ForceX];
            FyRFT{column, 5} = [FyRFT{column, 5} rft_result.ForceY];
            FzRFT{column, 5} = [FzRFT{column, 5} rft_result.ForceZ];
            
            depth{column, 5} = [depth{column, 5} -exp_result.avg_Z];
            depthRFT{column, 5} = [depthRFT{column, 5} -rft_result.depth];
            
        case 75
            
            Fx{column, 6} = [Fx{column, 6} -exp_result.avg_Fy];
            Fy{column, 6} = [Fy{column, 6} -exp_result.avg_Fx];
            Fz{column, 6} = [Fz{column, 6} -exp_result.avg_Fz];
            
            FxRFT{column, 6} = [FxRFT{column, 6} rft_result.ForceX];
            FyRFT{column, 6} = [FyRFT{column, 6} rft_result.ForceY];
            FzRFT{column, 6} = [FzRFT{column, 6} rft_result.ForceZ];
            
            depth{column, 6} = [depth{column, 6} -exp_result.avg_Z];
            depthRFT{column, 6} = [depthRFT{column, 6} -rft_result.depth];
            
        case 90
            
            Fx{column, 7} = [Fx{column, 7} -exp_result.avg_Fy];
            Fy{column, 7} = [Fy{column, 7} -exp_result.avg_Fx];
            Fz{column, 7} = [Fz{column, 7} -exp_result.avg_Fz];
            
            FxRFT{column, 7} = [FxRFT{column, 7} rft_result.ForceX];
            FyRFT{column, 7} = [FyRFT{column, 7} rft_result.ForceY];
            FzRFT{column, 7} = [FzRFT{column, 7} rft_result.ForceZ];
            
            depth{column, 7} = [depth{column, 7} -exp_result.avg_Z];
            depthRFT{column, 7} = [depthRFT{column, 7} -rft_result.depth];
        otherwise
            color = 'k';
    end
    
end


Fx_avg = zeros(9, 7);
Fy_avg = zeros(9, 7);
Fz_avg = zeros(9, 7);

Fx_err = zeros(9, 7);
Fy_err = zeros(9, 7);
Fz_err = zeros(9, 7);

Fx_rft_avg = zeros(9, 7);
Fy_rft_avg = zeros(9, 7);
Fz_rft_avg = zeros(9, 7);

depth_avg = zeros(9, 7);
depth_rft_avg = zeros(9, 7);

for i = 1:9
    for j = 1:7
        Fx_avg(i,j) = mean(Fx{i, j});
        Fy_avg(i,j) = mean(Fy{i, j});
        Fz_avg(i,j) = mean(Fz{i, j});
        
        Fx_err(i,j) = std(Fx{i, j})/sqrt(3);
        Fy_err(i,j) = std(Fy{i, j})/sqrt(3);
        Fz_err(i,j) = std(Fz{i, j})/sqrt(3);
        
        Fx_rft_avg(i,j) = mean(FxRFT{i, j});
        Fy_rft_avg(i,j) = mean(FyRFT{i, j});
        Fz_rft_avg(i,j) = mean(FzRFT{i, j});
        
        depth_avg(i,j) = mean(depth{i, j});
        depth_rft_avg(i,j) = mean(depthRFT{i, j});
    end
end

Fx_avg = Fx_avg';
Fy_avg = Fy_avg';
Fz_avg = Fz_avg';
Fx_err = Fx_err';
Fy_err = Fy_err';
Fz_err = Fz_err';

Fx_rft_avg = Fx_rft_avg';
Fy_rft_avg = Fy_rft_avg';
Fz_rft_avg = Fz_rft_avg';

depth_avg = depth_avg';
depth_rft_avg = depth_rft_avg';

deltaX = 0;
deltaY = 0;
deltaZ = 0;

R_xlist =[];
R_ylist =[];
R_zlist =[];


slip_angle = 5;
slip_ratio = 9;

err_slip = zeros(3, slip_angle);

for i = 1:slip_angle
    

    Xe = 0;
    Xm = 0;
    XeXm = 0;
    Xe2 = 0;
    Xm2 = 0;

    Ye = 0;
    Ym = 0;
    YeYm = 0;
    Ye2 = 0;
    Ym2 = 0;

    Ze = 0;
    Zm = 0;
    ZeZm = 0;
    Ze2 = 0;
    Zm2 = 0;
    for j = 1:slip_ratio
        
        exp_x = Fx_avg(i, j);
        exp_y = Fy_avg(i, j);
        exp_z = Fz_avg(i, j);
        
        rft_x = Fx_rft_avg(i,j);
        rft_y = Fy_rft_avg(i,j);
        rft_z = Fz_rft_avg(i,j);
        
        err_slip(1, i) = err_slip(1, i) + abs(exp_x - rft_x);
        err_slip(2, i) = err_slip(2, i) + abs(exp_y - rft_y);
        err_slip(3, i) = err_slip(3, i) + abs(exp_z - rft_z);
        
        deltaX = deltaX + abs(exp_x - rft_x);
        deltaY = deltaY + abs(exp_y - rft_y);
        deltaZ = deltaZ + abs(exp_z - rft_z);
        
        Xe = Xe + exp_x;
        Ye = Ye + exp_y;
        Ze = Ze + exp_z;

        Xm = Xm + rft_x;
        Ym = Ym + rft_y;
        Zm = Zm + rft_z;

        XeXm = XeXm + exp_x * rft_x;
        YeYm = YeYm + exp_y * rft_y;
        ZeZm = ZeZm + exp_z * rft_z;

        Xe2 = Xe2 + exp_x * exp_x;
        Ye2 = Ye2 + exp_y * exp_y;
        Ze2 = Ze2 + exp_z * exp_z;

        Xm2 = Xm2 + rft_x * rft_x;
        Ym2 = Ym2 + rft_y * rft_y;
        Zm2 = Zm2 + rft_z * rft_z;
    end
    err_slip(1, i) = err_slip(1, i)/slip_ratio;
    err_slip(2, i) = err_slip(2, i)/slip_ratio;
    err_slip(3, i) = err_slip(3, i)/slip_ratio;
    
    XR = (9 * XeXm - Xe * Xm)/(sqrt(9 * Xe2 - Xe * Xe) * sqrt(9 * Xm2 - Xm * Xm));
    YR = (9 * YeYm - Ye * Ym)/(sqrt(9 * Ye2 - Ye * Ye) * sqrt(9 * Ym2 - Ym * Ym));
    ZR = (9 * ZeZm - Ze * Zm)/(sqrt(9 * Ze2 - Ze * Ze) * sqrt(9 * Zm2 - Zm * Zm));
    R_xlist = [R_xlist XR];
    R_ylist = [R_ylist YR];
    R_zlist = [R_zlist ZR];
    
end

deltaX =deltaX / (slip_angle * slip_ratio)
deltaY =deltaY / (slip_angle * slip_ratio)
deltaZ =deltaZ / (slip_angle * slip_ratio)
R_xlist
R_ylist
R_zlist

figure()
title('depth')
hold on
for i = 1:slip_angle
    switch i
        case 1
            color = cmuColor('red-web');
        case 2
            color = cmuColor('gold');
        case 3
            color = cmuColor('teal');
        case 4
            color = cmuColor('sky-blue');
        case 5
            color = cmuColor('palladian-green');
        case 6
            color = cmuColor('blue-thread');
        case 7
            color = cmuColor('scots-rose');
        otherwise
            color = 'k';
    end
    
    
    e = plot(x, depth_avg(i,:),'o', 'MarkerEdgeColor', color);
    hold on
    p = plot(x, depth_rft_avg(i,:), 'Color', color);
end


figure()
sgtitle ('Wheel Forces Comparison between 3D-RFT+DEM and Experiment Results')
subplot(1,3,2)
title('Fy (Sidewall)')
xlim([-1.1 1])
hold on
xlabel('Slip Ratio')
ylabel('Force (N)')

subplot(1,3,1)
title('Fx (Tractive)')
xlim([-1.1 1])
hold on
xlabel('Slip Ratio')
ylabel('Force (N)')

subplot(1,3,3)
title('Fz (Load)')
xlim([-1.1 1])
hold on
xlabel('Slip Ratio')
ylabel('Force (N)')
for i = 1:slip_angle
    switch i
        case 1
            color = cmuColor('red-web');
        case 2
            color = cmuColor('gold');
        case 3
            color = cmuColor('teal');
        case 4
            color = cmuColor('sky-blue');
        case 5
            color = cmuColor('palladian-green');
        case 6
            color = cmuColor('blue-thread');
        case 7
            color = cmuColor('scots-rose');
        otherwise
            color = 'k';
    end
    subplot(1,3,1)
    e = errorbar(x,Fx_avg(i,:),Fx_err(i,:),'o');
    
    alpha = 0.2;
    e.Marker = 's';
    e.MarkerSize = 5;
    e.Color = color;
    e.CapSize = 10;
    e.LineWidth = 2;
    set([e.Bar, e.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e.Line.ColorData(1:3); 255*alpha])
    set(e.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [e.Cap.EdgeColorData(1:3); 255*alpha])
    
    p = plot(x, Fx_rft_avg(i,:), 'Color', color);
    
    subplot(1,3,2)
    e = errorbar(x,Fy_avg(i,:),Fy_err(i,:),'o');
    e.Marker = 's';
    e.MarkerSize = 5;
    e.Color = color;
    e.CapSize = 10;
    e.LineWidth = 2;
    set([e.Bar, e.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e.Line.ColorData(1:3); 255*alpha])
    set(e.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [e.Cap.EdgeColorData(1:3); 255*alpha])
    
    p = plot(x, Fy_rft_avg(i,:), 'Color', color);
    
    subplot(1,3,3)
    e = errorbar(x,Fz_avg(i,:),Fz_err(i,:),'o');
    
    e.Marker = 's';
    e.MarkerSize = 5;
    e.Color = color;
    e.CapSize = 10;
    e.LineWidth = 2;
    set([e.Bar, e.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e.Line.ColorData(1:3); 255*alpha])
    set(e.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [e.Cap.EdgeColorData(1:3); 255*alpha])
    
    p = plot(x, Fz_rft_avg(i,:), 'Color', color);
end


% for i=1:length(RFToutput)
%     rft_result = RFToutput(i);
%     %Select color of plotted point
%     switch rft_result.beta
%         case 0
%             color = cmuColor('red-web');
%         case 15
%             color = cmuColor('gold');
%         case 30
%             color = cmuColor('teal');
%         case 45
%             color = cmuColor('sky-blue');
%         case 60
%             color = cmuColor('palladian-green');
%         case 75
%             color = cmuColor('blue-thread');
%         case 90
%             color = cmuColor('scots-rose');
%         otherwise
%             color = 'k';
%     end
%     
%     subplot(1,3,1)
%     p = plot(rft_result.slip, rft_result.ForceX, 's', 'MarkerEdgeColor', color);
%     p.Color(4) = 0.25;
%     
%     subplot(1,3,2)
%     p = plot(rft_result.slip, rft_result.ForceY, 's', 'MarkerEdgeColor', color);
%     p.Color(4) = 0.25;
%     
%     subplot(1,3,3)
%     p = plot(rft_result.slip, rft_result.ForceZ, 's', 'MarkerEdgeColor', color);
%     p.Color(4) = 0.25;
% end

leg(1) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('red-web'));
leg(2) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('gold'));
leg(3) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('teal'));
leg(4) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('sky-blue'));
leg(5) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('palladian-green'));
legg = legend(leg,'0', '15', '30', '45', '60');
title(legg, 'Slip Angle');
hold off

function [column] = switchstuff(wr)
    
    switch wr
        case 0
            column = 1;
        case 3
            column = 2;
        case 5
            column = 3;
        case 8
            column = 4;
        case 10
            column = 5;
        case 12
            column = 6;
        case 20
            column = 7;
        case 40
            column = 8;
        case 100
            column = 9;
        otherwise
            column = 0;
           
    end
end
