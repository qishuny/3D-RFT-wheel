%% gives you 4x4 HT, (for 2.5D) [x,y,z,1]'
% there is no pitch or roll. Only yaw
function HT = get_HT4(x,y,z,th)
HT = [cos(th) -sin(th) 0 x;
    sin(th) cos(th) 0 y;
    0 0 1 z;
    0 0 0 1];
end