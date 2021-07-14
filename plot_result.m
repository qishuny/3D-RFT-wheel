close all
load('output/RFToutput.mat')

%% Forces Figure
figure()
sgtitle ('Wheel Forces')
subplot(1,3,2)
title('Fy (Sidewall)')
axis([-1.1 1 -150 50])
hold on
xlabel('Slip Ratio')
ylabel('Force (N)')

subplot(1,3,1)
title('Fx (Tractive)')
axis([-1.1 1 -150 175])
hold on
xlabel('Slip Ratio')
ylabel('Force (N)')

subplot(1,3,3)
title('Fz (Load)')
axis([-1.1 1 0 200])
hold on
xlabel('Slip Ratio')
ylabel('Force (N)')

for i=1:length(RFToutput)
    result = RFToutput(i);
    %Select color of plotted point
    switch result.beta
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
    
    subplot(1,3,2)
    plot(result.slip, result.ForceY, 'o', 'MarkerEdgeColor', color)

    subplot(1,3,1)
    plot(result.slip, result.ForceX, 'o', 'MarkerEdgeColor', color)
    
    subplot(1,3,3)
    plot(result.slip, result.ForceZ, 'o', 'MarkerEdgeColor', color)
    
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




