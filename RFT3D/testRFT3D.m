% load experiment data
load('../output/all_smooth_data_2.mat')

% load wheel point data
wheeldata = matfile('../data/smooth_wheel_125.mat');
radius = 62.5;
vcenter = 10;
scale = 0.5;

h = waitbar(0,'Initializing waitbar...');
for i=1:length(all_results)
    result = all_results(i);
    wr = result.Vry;
    w = wr / radius;
    sinkage = abs(result.avg_Z);
    slipAngle = result.beta * pi / 180;
    

    [Force] = RFT3Dfunc(wheeldata, radius, slipAngle, w, vcenter, sinkage, scale);
    Fsidewall = Force(1);
    Ftractive = -Force(2);
    Fload = Force(3);
    RFToutput(i) = struct('ForceX', Ftractive, ...
        'ForceY', Fsidewall , ...
        'ForceZ', Fload, ...
        'wr', result.Vry, ...
        'depth', result.avg_Z, ...
        'beta', result.beta, ...
        'slip', result.slip); 
    
    waitbar(i/length(all_results), h, 'In progress...')
end
waitbar(1,h,'Completed.');
disp("Done.");

close(h);
save('../output/RFTtestoutput.mat','RFToutput');


