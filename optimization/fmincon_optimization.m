%Run tuning_optimization_wrapper and call this function to tune the
%values for a and b on the data. 
function [x, fval, exitflag] = fmincon_optimization(all_results)

% Optimization-based tuning of sf1 sf2 

%Constraints and bounds
A = [];
b = []'; % Ax <= b -> a0 + a1 <=1, -b0 - b1 <= 1
lb = [0,0];
ub = [1,1];
Aeq = []; % 
beq = []'; % Ax = b -> a0*theta_f0 - theta_m = 0, b0*theta_f0 - theta_r = 0
nonlcon = []; %No nonlinear constraints

x0 = [0.9 0.9]'; %Starting guesses, x = [a0 a1 b0 b1 theta_m theta_r n k K c]'
%options = optimoptions(@fmincon, 'Algorithm', 'sqp', 'Display', 'iter','OptimalityTolerance', 1e-16, 'FunctionTolerance', 1e-16, 'StepTolerance', 1e-16, 'DiffMinChange', 1e-16, 'FiniteDifferenceStepSize', 1e-6); %Saved as a comment as a record of the settings I have tried playing with
options = optimoptions(@fmincon, 'Algorithm', 'interior-point', 'Display', 'iter', 'ConstraintTolerance', 1e-12, 'OptimalityTolerance', 1e-3, 'StepTolerance', 1e-16); %Sets the options for fmincon, algorithm choice shouldn't affect much
[x, fval, exitflag] = fmincon(@(x)tuning_objective_fn(x, all_results), x0, A, b, Aeq, beq, lb, ub, nonlcon, options);

