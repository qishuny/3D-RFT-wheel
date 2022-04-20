function H_new = mtx_method_update_height(H,dx,blade_idx,dt)
    a = 1/dx * (dx^2 * 1/2^3); %0.001 / dx * dt; %0.1*0.001 / dx; %sand flow rate
    b = tand(29); %0.6; %threshold slope
    %%
    if nargin > 2 %when there is blade
        
    % flow from right
    H_right = [H(:,2:end) H(:,end)];
%     for i = 1:size(blade_idx,2) % blocking flow from blade locations
%         
%         H_right(blade_idx(1,i),blade_idx(2,i)) = H(blade_idx(1,i),blade_idx(2,i)); %index at blade
%         try
%         H_right(blade_idx(1,i),blade_idx(2,i)-1) = H(blade_idx(1,i),blade_idx(2,i)-1); %index adjacent to blade
%         catch
%         end
%     end
    slope = (H_right - H)/dx;
%     H_update1 = temp2.*(abs(temp2) > b)*alphadx;
%     H_update1 = sign(temp2).*(abs(temp2) - b)*alphadx;
    dH1 = get_update(slope,b,a); % (slope +- thres)*(abs(slope) > thres)
    
    % flow from left
    H_next = [H(:,1) H(:,1:end-1)];
%     for i = 1:size(blade_idx,2)
%         
%         H_next(blade_idx(1,i),blade_idx(2,i)) = H(blade_idx(1,i),blade_idx(2,i));%index at blade
%         try
%         H_next(blade_idx(1,i),blade_idx(2,i)+1) = H(blade_idx(1,i),blade_idx(2,i)+1);%index adjacent to blade
%         catch
%         end
%     end   
    slope = (H_next - H)/dx;
%     H_update2 = temp2.*(abs(temp2) > b)*alphadx;
%     H_update2 = sign(temp2).*(abs(temp2) - b)*alphadx;
    dH2 = get_update(slope,b,a);

    % flow from up
    H_next = [H(1,:); H(1:end-1,:)];
%     for i = 1:size(blade_idx,2)
%         
%         H_next(blade_idx(1,i),blade_idx(2,i)) = H(blade_idx(1,i),blade_idx(2,i));%index at blade
%         try
%         H_next(blade_idx(1,i)+1,blade_idx(2,i)) = H(blade_idx(1,i)+1,blade_idx(2,i));%index adjacent to blade
%         catch
%         end
%     end 
    slope = (H_next - H)/dx;
    dH3 = get_update(slope,b,a);

    % flow from down
    H_next = [H(2:end,:) ; H(end,:)];
%     for i = 1:size(blade_idx,2)
% 
%         H_next(blade_idx(1,i),blade_idx(2,i)) = H(blade_idx(1,i),blade_idx(2,i));
%         try
%         H_next(blade_idx(1,i)-1,blade_idx(2,i)) = H(blade_idx(1,i)-1,blade_idx(2,i));
%         catch
%         end
%     end 
    slope = (H_next - H)/dx;
    dH4 = get_update(slope,b,a);

%     H_new = H + dH1 + dH2 + dH3 + dH4;
    
    side_area = 0.431884367759256;
    corner_area = 0.242402684527083;
    r = corner_area/side_area; %the ratio between side and corner area

    % flow from right-up
    H_next = [H(1,1:end-1);H(1:end-1,2:end)];
    H_next = [H_next H(:,end)];
%     for i = 1:size(blade_idx,2)
% 
%         H_next(blade_idx(1,i),blade_idx(2,i)) = H(blade_idx(1,i),blade_idx(2,i));
%         try
%         H_next(blade_idx(1,i)+1,blade_idx(2,i)-1) = H(blade_idx(1,i)+1,blade_idx(2,i)-1);
%         catch
%         end
%     end     
    slope = (H_next - H)/(sqrt(2)*dx);
    dH5 = get_update(slope,b,a)*r;
    
    % flow from right-down
    H_next = [H(2:end,2:end) H(1:end-1,end)];
    H_next = [H_next; H(end,:)];
%     for i = 1:size(blade_idx,2)
% 
%         H_next(blade_idx(1,i),blade_idx(2,i)) = H(blade_idx(1,i),blade_idx(2,i));
%         try
%         H_next(blade_idx(1,i)-1,blade_idx(2,i)-1) = H(blade_idx(1,i)-1,blade_idx(2,i)-1);
%         catch
%         end
%     end 
    slope = (H_next - H)/(sqrt(2)*dx);
    dH6 = get_update(slope,b,a)*r;
    
    % flow from left-up
    H_next = [H(1,2:end); H(1:end-1,1:end-1)];
    H_next = [H(:,1) H_next];
%     for i = 1:size(blade_idx,2)
% 
%         H_next(blade_idx(1,i),blade_idx(2,i)) = H(blade_idx(1,i),blade_idx(2,i));
%         try
%         H_next(blade_idx(1,i)+1,blade_idx(2,i)+1) = H(blade_idx(1,i)+1,blade_idx(2,i)+1);
%         catch
%         end
%     end 
    slope = (H_next - H)/(sqrt(2)*dx);
    dH7 = get_update(slope,b,a)*r;
    
    % flow from left-down
    H_next = [H(1:end-1,1) H(2:end,1:end-1)];
    H_next = [H_next; H(end,:)];
%     for i = 1:size(blade_idx,2)
% 
%         H_next(blade_idx(1,i),blade_idx(2,i)) = H(blade_idx(1,i),blade_idx(2,i));
%         try
%         H_next(blade_idx(1,i)-1,blade_idx(2,i)+1) = H(blade_idx(1,i)-1,blade_idx(2,i)+1);
%         catch
%         end
%     end 
    slope = (H_next - H)/(sqrt(2)*dx);
    dH8 = get_update(slope,b,a)*r;
       
    
          
    no_flow_idx = dH1(blade_idx) < 0;
    dH1(no_flow_idx) = 0;
    no_flow_idx = dH2(blade_idx) < 0;
    dH2(no_flow_idx) = 0;
    no_flow_idx = dH3(blade_idx) < 0;
    dH3(no_flow_idx) = 0;
    no_flow_idx = dH4(blade_idx) < 0;
    dH4(no_flow_idx) = 0;
    no_flow_idx = dH5(blade_idx) < 0;
    dH5(no_flow_idx) = 0;
    no_flow_idx = dH6(blade_idx) < 0;
    dH6(no_flow_idx) = 0;
    no_flow_idx = dH7(blade_idx) < 0;
    dH7(no_flow_idx) = 0;
    no_flow_idx = dH8(blade_idx) < 0;
    dH8(no_flow_idx) = 0;
    H_new = H + dH1 + dH2 + dH3 + dH4 ...
              + dH5 + dH6 + dH7 + dH8;
     
    elseif nargin == 2 % when there is no blade
     
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
    dH7 = get_update(slope,b,a)*r;get_update
    
    % flow from left-down
    H_next = [H(1:end-1,1) H(2:end,1:end-1)];
    H_next = [H_next; H(end,:)];
    slope = (H_next - H)/(sqrt(2)*dx);
%     dH8 = temp2.*(abs(temp2)>b)*alphadx*r;
    dH8 = get_update(slope,b,a)*r;
    
    H_new = H + dH1 + dH2 + dH3 + dH4 ...
              + dH5 + dH6 + dH7 + dH8;  
    
    end 
    
    function update = get_update(slope,thres,a)
%         update = slope.*(abs(slope)>thres)*adx;
        update = sign(slope).*(abs(slope) - thres).*(abs(slope)>thres)*a;
    end

end