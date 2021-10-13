load('output/all_smooth_data_2.mat')
% load('output/RFToutput.mat')
load('output/RFTDEMoutput.mat')


n = length(all_results);

deltaX = 0;
deltaY = 0;
deltaZ = 0;

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

for i=1:n
    exp_result = all_results(i);
    rft_result = RFToutput(i);
    
    
    exp_x = -exp_result.avg_Fy;
    exp_y = -exp_result.avg_Fx;
    exp_z = -exp_result.avg_Fz;
    
    rft_x = rft_result.ForceX;
    rft_y = rft_result.ForceY;
    rft_z = rft_result.ForceZ;
    
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

deltaX = deltaX / length(all_results);
deltaY = deltaY / length(all_results);
deltaZ = deltaZ / length(all_results);

%%
XR = (n * XeXm - Xe * Xm)/(sqrt(n * Xe2 - Xm * Xm) * sqrt(n * Xm2 - Xe * Xe))
YR = (n * YeYm - Ye * Ym)/(sqrt(n * Ye2 - Ym * Ym) * sqrt(n * Ym2 - Ye * Ye))
ZR = (n * ZeZm - Ze * Zm)/(sqrt((n * Ze2) - (Zm * Zm)) * sqrt((n * Zm2) - (Ze * Ze)))


