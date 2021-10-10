H = [ 0.59591124,  0.54325211,  0.54608145;
      0.63882917, -0.26234925, -0.82771596;
      0.78473912,  0.40098575,  0.0470053];
dx = 0.01;
a = 1/dx * dx^2 * 1/2^3; %0.001 / dx * dt; %0.1*0.001 / dx; %sand flow rate
b = tand(29); %0.6; %threshold slope
%%
% flow from right
H_next = [H(:,2:end) H(:,end)];
slope = (H_next - H)/dx;
%     dH1 = temp2.*(abs(temp2) > b)*alphadx;
dH1 = get_update(slope,b,a);

% flow from left
H_next = [H(:,1) H(:,1:end-1)];
slope = (H_next - H)/dx;
%     dH2 = temp2.*(abs(temp2) > b)*alphadx;
dH2 = get_update(slope,b,a);

% flow from up
H_next = [H(1,:); H(1:end-1,:)];
slope = (H_next - H)/dx;
%     dH3 = temp2.*(abs(temp2) > b)*alphadx;
dH3 = get_update(slope,b,a);

% flow from down
H_next = [H(2:end,:) ; H(end,:)];
slope = (H_next - H)/dx;
%     dH4 = temp2.*(abs(temp2) > b)*alphadx;
dH4 = get_update(slope,b,a);

%     H_new = H + dH1 + dH2 + dH3 + dH4;

side_area = 0.431884367759256;
corner_area = 0.242402684527083;
r = corner_area/side_area; %the ratio between side and corner area

% flow from right-up
H_next = [H(1,1:end-1);H(1:end-1,2:end)];
H_next = [H_next H(:,end)];
slope = (H_next - H)/(sqrt(2)*dx);
%     dH5 = temp2.*(abs(temp2)>b)*alphadx*r;
dH5 = get_update(slope,b,a)*r;

% flow from right-down
H_next = [H(2:end,2:end) H(1:end-1,end)];
H_next = [H_next; H(end,:)];
slope = (H_next - H)/(sqrt(2)*dx);
%     dH6 = temp2.*(abs(temp2)>b)*alphadx*r;
dH6 = get_update(slope,b,a)*r;

% flow from left-up
H_next = [H(1,2:end); H(1:end-1,1:end-1)];
H_next = [H(:,1) H_next];
slope = (H_next - H)/(sqrt(2)*dx);
%     dH7 = temp2.*(abs(temp2)>b)*alphadx*r;
dH7 = get_update(slope,b,a)*r;

% flow from left-down
H_next = [H(1:end-1,1) H(2:end,1:end-1)];
H_next = [H_next; H(end,:)];
slope = (H_next - H)/(sqrt(2)*dx);
%     dH8 = temp2.*(abs(temp2)>b)*alphadx*r;
dH8 = get_update(slope,b,a)*r;

H_new = H + dH1 + dH2 + dH3 + dH4 ...
    + dH5 + dH6 + dH7 + dH8;
    

function update = get_update(slope,thres,a)
%         update = slope.*(abs(slope)>thres)*adx;
update = sign(slope).*(abs(slope) - thres).*(abs(slope)>thres)*a;
end