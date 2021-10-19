%% load data

load('/home/qishuny/3D-RFT-wheel/output/forceBalance.mat')
load('/home/qishuny/3D-RFT-wheel/output/all_smooth_data_2.mat')

bList = [0 15 30 45 60 75 90];
wlist = [0 3 5 8 10 12 20 40 100];
slipList = [-1 -0.7 -0.5 -0.2 0 0.23 0.5 0.75 0.89];
%%

rft_Fx = zeros(1, length(all_results));
rft_Fy = zeros(1, length(all_results));
rft_Fz = zeros(1, length(all_results));

count = 0;
for i = 1:length(all_results)
    
    exp_result = all_results(i);
    wr_exp = exp_result.Vry;
    beta = exp_result.beta;
    idxi = find(bList == beta)
    idxj = find(wlist == wr_exp)
    targetFz = -exp_result.avg_Fz;
    targetFx = -exp_result.avg_Fy;
    targetFy = -exp_result.avg_Fx;
    Fz = 0;
    Fx = 0;
    Fy = 0;
    for j = 1:1000
        
        if abs(Depth(idxj, idxi, j).ForceZ - targetFz) < 0.1
            Fz =Depth(idxj, idxi, j).ForceZ;
            Fx =Depth(idxj, idxi, j).ForceX;
            Fy =Depth(idxj, idxi, j).ForceY;
            count = count + 1;
            rft_Fx(i) = Fx;
            rft_Fy(i) = Fy;
            rft_Fz(i) = Fz;
        end
    end
    RFToutput(i) = struct('ForceX', Fx, ...
        'ForceY', Fy , ...
        'ForceZ', Fz, ...
        'wr', wr_exp, ...
        'depth', exp_result.avg_Z, ...
        'beta', beta, ...
        'slip', exp_result.slip);

end
%%


save('output/RFTbalanceoutput.mat','RFToutput');