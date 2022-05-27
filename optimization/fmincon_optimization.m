%Run tuning_optimization_wrapper and call this function to tune the
%values for a and b on the data. 
function [x, fval, exitflag] = fmincon_optimization(wheeldata, all_results)

% Optimization-based tuning of sf1 sf2 

%Constraints and bounds
A = [];
b = []'; % Ax <= b -> a0 + a1 <=1, -b0 - b1 <= 1
lb = [0];
ub = [1];
Aeq = []; % 
beq = []'; % Ax = b -> a0*theta_f0 - theta_m = 0, b0*theta_f0 - theta_r = 0
nonlcon = []; %No nonlinear constraints

% x0 = [0.5 0.5]'; %Starting guesses, x = [a0 a1 b0 b1 theta_m theta_r n k K c]'
x0 = [0.3]';
options = optimoptions(@fmincon, 'Algorithm', 'interior-point', 'Display', 'iter', 'ConstraintTolerance', 1e-12, 'OptimalityTolerance', 1e-3, 'StepTolerance', 1e-16); %Sets the options for fmincon, algorithm choice shouldn't affect much
[x, fval, exitflag] = fmincon(@(x)tuning_objective_fn(wheeldata, x, all_results), x0, A, b, Aeq, beq, lb, ub, nonlcon, options);


