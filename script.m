%% load data

load('/home/qishuny/3D-RFT-wheel/output/forceBalance.mat')
load('/home/qishuny/3D-RFT-wheel/output/all_smooth_data_2.mat')

bList = [0 15 30 45 60 75 90];
wlist = [0 3 5 8 10 12 20 40 100];
slipList = [-1 -0.7 -0.5 -0.2 0 0.23 0.5 0.75 0.89];


%%
for i = 1:length(all_results)
    exp_result = all_results(1);
    wr_exp = exp_result.Vry;
    beta = exp_result.beta;
    idxi = find(bList == beta);
    idxj= find(wlist == wr_exp);
end
