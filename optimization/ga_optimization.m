%Run tuning_optimization_wrapper and call this function to tune the
%values for a and b on the data. 
function [X, Y, Z, x] = ga_optimization(all_results)

% Optimization-based tuning of sf1 sf2 
sf1 = linspace(0.1,0.99,5);
sf2 = linspace(0.1,0.99,5);
[X,Y] = meshgrid(sf1, sf2);

n = size(X,1) * size(X,2);


for i = 1:n
    
    xtemp = [X(i) Y(i)];
    Z(i) = tuning_objective_fn(xtemp, all_results);
    
end
Z = reshape(Z,size(X))

surf(X,Y,Z,'MeshStyle','none')
colormap 'jet'
view(-26,43)
xlabel('x(1)')
ylabel('x(2)')
title('ps\_example(x)')


%Constraints and bounds
A = [];
b = []'; % Ax <= b -> a0 + a1 <=1, -b0 - b1 <= 1
lb = [0,0];
ub = [1,1];
Aeq = []; % 
beq = []'; % Ax = b -> a0*theta_f0 - theta_m = 0, b0*theta_f0 - theta_r = 0
nonlcon = []; %No nonlinear constraints

x = ga(fun,2,A,b)

