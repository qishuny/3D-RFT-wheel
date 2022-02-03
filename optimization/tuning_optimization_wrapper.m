load('output/all_smooth_data_2.mat');
wheeldata = matfile('data/smooth_wheel_125.mat');

[x, fval, exitflag] = fmincon_optimization(wheeldata, all_results);



