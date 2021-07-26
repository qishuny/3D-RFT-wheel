%%
function indices = get_coord(p,dx,minx,miny)
if nargin == 2
    x_idx = round(p(1,:)/dx)+1;
    y_idx = round(p(2,:)/dx)+1;
    indices = [x_idx; y_idx];
elseif nargin ==4
    x_vec = p(1,:) - minx;
    y_vec = p(2,:) - miny;
    x_idx = round(x_vec/dx)+1;
    y_idx = round(y_vec/dx)+1;
    indices = [x_idx; y_idx];
end