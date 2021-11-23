load('output/all_smooth_data_2.mat')
wheeldata = matfile('data/smooth_wheel_125.mat');


% optimized for sf1 sf2
% sf1 = 0.175;
% sf2 = 0.5;

%optimized for sf1
sf1 = 0.075;
sf2 = 0.6975;

h = waitbar(0,'Initializing waitbar...');
for i=1:length(all_results)
    result = all_results(i);
    wr = result.Vry;
    if wr == 0
        wr = 0.00001;
    end  
    sinkage = abs(result.avg_Z);

    slipAngle = result.beta * pi / 180;
    

    [Ftractive, Fsidewall, Fload] = RFT3DDEMfunc(wheeldata, slipAngle, wr, sinkage, sf1, sf2);

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
save('output/RFTDEMoptoutput.mat','RFToutput');