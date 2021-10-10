function update = get_update___(slope,thres,a)
%         update = slope.*(abs(slope)>thres)*adx;
        update = sign(slope).*(abs(slope) - thres).*(abs(slope)>thres)*a;
end