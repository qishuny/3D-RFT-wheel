load('../output/all_smooth_data_2.mat');
wheeldata = matfile('../data/smooth_wheel_125.mat');

x = linspace(0, 1, 10);
for i = 1:10
    i
    error(i) = tuning_objective_fn(wheeldata, x(i), all_results);
end

plot(x, error)