%Run tuning_optimization_wrapper and call this function to tune the
%values for a and b on the data. 
function [X, Y, Z, x] = ga_optimization(wheeldata, all_results)

% Optimization-based tuning of sf1 sf2 

% sf1 = linspace(0.01,0.5,5);
% sf2 = linspace(0.6975,0.6975,5);
% [X,Y] = meshgrid(sf1, sf2);
% 
% n = size(X,1) * size(X,2);
% 
% 
% for i = 1:n
%     
%     xtemp = [X(i) Y(i)];
%     Z(i) = tuning_objective_fn(wheeldata, xtemp, all_results);
%     i
% end
% Z = reshape(Z,size(X))
% 
% surf(X,Y,Z,'MeshStyle','none')
% colormap 'jet'
% view(-26,43)
% xlabel('x(1)')
% ylabel('x(2)')
% title('ps\_example(x)')

% tune sf1

sf1 = linspace(0.01,0.99,5);
sf2 = 0.6975;
for i = 1:5
    xtemp = [sf1(i) sf2];
    Z(i) = tuning_objective_fn(wheeldata, xtemp, all_results);
end
plot(sf1, Z)
X = sf1;
Y = sf2;

%Constraints and bounds
A = [];
b = []'; % Ax <= b -> a0 + a1 <=1, -b0 - b1 <= 1
lb = [0,0];
ub = [1,1];
Aeq = []; % 
beq = []'; % Ax = b -> a0*theta_f0 - theta_m = 0, b0*theta_f0 - theta_r = 0
nonlcon = []; %No nonlinear constraints
x = 0;
% x = ga(fun,2,A,b)

