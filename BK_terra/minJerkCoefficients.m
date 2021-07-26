%%
function a = minJerkCoefficients(x0,x0d,x0dd,xf,xfd,xfdd,tf)
% get coefficients of minJerk trajectory

% a0 = x0;
% a1 = x0d;
% a2 = 0.5*x0dd;
% 
% tfMatrix = [ tf^3 tf^4 tf^5; 3*tf^2 4*tf^3 5*tf^4; 6*tf 12*tf^2 20*tf^3];
% BC = [xf-a0-a1-a2; xfd-a1-2*a2; xfdd-2*a2];
% a345 = tfMatrix*BC;
% 
% a = [a0 a1 a2 a345(1) a345(2) a345(3)];
T = tf;
T2 = T*T; T3 = T2*T;
T4 = T3*T; T5= T4*T;
a = zeros(6,1);
a(1) = x0;
a(2) = x0d;
a(3) = x0dd/2;
b= [T3 T4 T5 ; 3*T2 4*T3 5*T4; 6*T 12*T2 20*T3];
c = [ xf - a(1) - a(2)*T - a(3)*T2; xfd - a(2) - 2*a(3)*T;
xfdd - 2*a(3)];
a(4:6,1)=pinv(b)*c;
end