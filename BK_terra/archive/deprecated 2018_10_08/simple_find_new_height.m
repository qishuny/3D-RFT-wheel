function H_new = simple_find_new_height(H,dx,alpha,dt,angle_of_repose)
%     angle_of_repose = 0.6; %threshold slope
    %% 
    % flow from right
    H_right = [H(:,2:end) H(:,end)];
    slope = (H_right - H)/dx;

    H_update1 = get_update(slope,angle_of_repose,alpha,dx,dt);
    
    % flow from left
    H_next = [H(:,1) H(:,1:end-1)];
    slope = (H_next - H)/dx;

    H_update2 = get_update(slope,angle_of_repose,alpha,dx,dt);
    H_new = H + H_update1 + H_update2;

    function update = get_update(slope,thres,alpha,dx,dt)
        update = alpha*dx*sign(slope).*(abs(slope) - thres).*(abs(slope)>thres)*dt/dx^2;
    end
end