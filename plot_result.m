close all
load('output/RFToutput.mat')
load('output/all_smooth_data_2.mat')
%% RFT output Figure
figure()
sgtitle ('RFT Wheel Forces')
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

for i=1:length(RFToutput)
    rft_result = RFToutput(i);
    %Select color of plotted point
    switch rft_result.beta
        case 0
            color = cmuColor('red-web');
        case 15
            color = cmuColor('gold');
        case 30
            color = cmuColor('teal');
        case 45
            color = cmuColor('sky-blue');
        case 60
            color = cmuColor('palladian-green');
        case 75
            color = cmuColor('blue-thread');
        case 90
            color = cmuColor('scots-rose');
        otherwise
            color = 'k';
    end
    
    subplot(1,3,1)
    plot(rft_result.slip, rft_result.ForceX, 'o', 'MarkerEdgeColor', color)
    
    subplot(1,3,2)
    plot(rft_result.slip, rft_result.ForceY, 'o', 'MarkerEdgeColor', color)
    
    subplot(1,3,3)
    plot(rft_result.slip, -rft_result.ForceZ, 'o', 'MarkerEdgeColor', color)
    
end
leg(1) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('red-web'));
leg(2) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('gold'));
leg(3) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('teal'));
leg(4) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('sky-blue'));
leg(5) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('palladian-green'));
leg(6) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('blue-thread'));
leg(7) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('scots-rose')); 
legg = legend(leg,'0', '15', '30', '45', '60', '75', '90');
title(legg, 'Slip Angle');
hold off

%% Experiment Figure
figure()
sgtitle ('EXP Wheel Forces')
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

for i=1:length(all_results)
    exp_result = all_results(i);
    %Select color of plotted point
    switch exp_result.beta
        case 0
            color = cmuColor('red-web');
        case 15
            color = cmuColor('gold');
        case 30
            color = cmuColor('teal');
        case 45
            color = cmuColor('sky-blue');
        case 60
            color = cmuColor('palladian-green');
        case 75
            color = cmuColor('blue-thread');
        case 90
            color = cmuColor('scots-rose');
        otherwise
            color = 'k';
    end
    
    subplot(1,3,1)
    plot(exp_result.slip, -exp_result.avg_Fy, 'o', 'MarkerEdgeColor', color)
    
    subplot(1,3,2)
    plot(exp_result.slip, -exp_result.avg_Fx, 'o', 'MarkerEdgeColor', color)
    
    subplot(1,3,3)
    plot(exp_result.slip, exp_result.avg_Fz, 'o', 'MarkerEdgeColor', color)
    
end
leg(1) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('red-web'));
leg(2) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('gold'));
leg(3) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('teal'));
leg(4) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('sky-blue'));
leg(5) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('palladian-green'));
leg(6) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('blue-thread'));
leg(7) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('scots-rose')); 
legg = legend(leg,'0', '15', '30', '45', '60', '75', '90');
title(legg, 'Slip Angle');
hold off

%% plot percentage error
slipList = zeros(size(all_results,1), 1);
angleList = zeros(size(all_results,1), 1);

exp_FxList = zeros(size(all_results,1), 1);
exp_FyList = zeros(size(all_results,1), 1);
exp_FzList = zeros(size(all_results,1), 1);

rft_FxList = zeros(size(all_results,1), 1);
rft_FyList = zeros(size(all_results,1), 1);
rft_FzList = zeros(size(all_results,1), 1);

per_errX = zeros(size(all_results,1), 1);
per_errY = zeros(size(all_results,1), 1);
per_errZ = zeros(size(all_results,1), 1);

errX = zeros(size(all_results,1), 1);
errY = zeros(size(all_results,1), 1);
errZ = zeros(size(all_results,1), 1);
for i = 1:size(all_results,1)
    slipList(i) = all_results(i).slip;
    angleList(i) = all_results(i).beta;
    
    exp_FxList(i) = -all_results(i).avg_Fy;
    exp_FyList(i) = -all_results(i).avg_Fx;
    exp_FzList(i) = -all_results(i).avg_Fz;
    
    rft_FxList(i) = RFToutput(i).ForceX;
    rft_FyList(i) = RFToutput(i).ForceY;
    rft_FzList(i) = RFToutput(i).ForceZ;
    
    errX(i) = rft_FxList(i) - exp_FxList(i);
    errY(i) = rft_FyList(i) - exp_FyList(i);
    errZ(i) = rft_FzList(i) - exp_FzList(i);
    
    
    per_errX(i) = (rft_FxList(i) - exp_FxList(i)) / exp_FxList(i) * 100;
    per_errY(i) = (rft_FyList(i) - exp_FyList(i)) / exp_FyList(i) * 100;
    per_errZ(i) = (rft_FzList(i) - exp_FzList(i)) / exp_FzList(i) * 100;
    
end

figure

subplot(1,3,1)

scatter3(slipList, angleList, errY)
subplot(1,3,2)
scatter3(slipList, angleList, errY)
subplot(1,3,3)
scatter3(slipList, angleList, errZ)
title('force error')

% figure
% subplot(1,3,1)
% scatter3(slipList, angleList, per_errX)
% subplot(1,3,2)
% scatter3(slipList, angleList, per_errY)
% subplot(1,3,3)
% scatter3(slipList, angleList, per_errZ)
% title('percentage error')

