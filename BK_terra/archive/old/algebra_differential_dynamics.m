%% Practicing Holonomic Constraint

m = 10; % kg
c = [1 0 3 0]; %[1 0 1 0]; %[0 0 1 0]; % gradient

V = eye(4) - c'*(c*c')^(-1)*c;
A = [0 1 0 0;
     0 0 0 0;
     0 0 0 1;
     0 0 0 0];
B = [0 1/m 0 0]';

u = 1; % N

x = zeros(4,10);
dT = 0.1; % s
T = zeros(1,10);
for i = 1:10
    dx = V*A*x(:,i) + V*B*u;
    x(:,i+1) = x(:,i) + dx*dT;
    T(i) = dT*(i-1);
end
plot(x(1,:),x(3,:),'b');





