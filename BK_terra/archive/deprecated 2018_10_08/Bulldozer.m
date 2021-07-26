%%
classdef Bulldozer < handle
    properties
        body_pos;
        body_ori;
        rel_pos_blade;
        rel_pos_blade_front;
        rel_pos_blade_back;
        blade_pos;
        blade_pos_front;
        blade_pos_back;
        dx;
        trans_matrix;
        
        x; y; th; %dx;
        Blade;
        BladeFront;
    end
    

    methods 
        
        function this = Bulldozer()
            this.dx = 0.1;
            this.x = 0; this.y = 0; this.th = 0; 
            this.Blade = ContactPart();
            this.BladeFront = ContactPart(0,0.2,0,this.dx);
            
            this.body_pos = [1; 0.3];
            % blade's x coord (first row)
            % blade's y coord (second row)
            this.body_ori = 0; % in rad
            th = this.body_ori;
            % transform matrix Rwb, bull frame to world  frame
            this.trans_matrix = [cos(th) -sin(th) this.body_pos(1);
                                 sin(th)  cos(th) this.body_pos(2);
                                 0        0       1               ];
            
            % relative blade pos
            % p_blade wrt bulldozer
            % sand flow from these positions are blocked
            this.rel_pos_blade = [-0.5 -0.4 -0.3 -0.2 -0.1 0   0.1 0.2 0.3 0.4 0.5;
                              0.2   0.2  0.2  0.2  0.2 0.2 0.2 0.2 0.2 0.2 0.2;
                              1     1    1    1    1   1   1   1   1   1   1  ];
                              % row of ones is for scaling factor (for
                              % tranform matrix
                              
            % relative pos of blade back
            % sand flow  from these positions are blocked as well
            this.rel_pos_blade_back = this.rel_pos_blade + ...
                                   [0;-this.dx;0]*ones(1,length(this.rel_pos_blade));                             
            % relative blade pos front info
            % front of the blade where sand will be displaced to             
            this.rel_pos_blade_front = this.rel_pos_blade + ...
                                   [0;this.dx;0]*ones(1,length(this.rel_pos_blade));
                               
%             this.rel_pos_blade_front = [-0.5 -0.4 -0.3 -0.2 -0.1 0   0.1 0.2 0.3 0.4 0.5;
%                                          0.4  0.4  0.4  0.4  0.4 0.4 0.4 0.4 0.4 0.4 0.4;
%                                          1     1    1    1    1   1   1   1   1   1   1  ];
            
            % absolute blade pos 
            this.blade_pos = this.trans_matrix*this.rel_pos_blade;
            this.blade_pos_front = this.trans_matrix*this.rel_pos_blade_front;
            this.blade_pos_back = this.trans_matrix*this.rel_pos_blade_back;                  
                 
        end

        function synchronize(this,dx,pos_x,pos_y)
            % synchronize with Heightmap 
            % have the same dx (so that my blade info match
            this.dx = dx;
            this.body_pos = [pos_x; pos_y];
            this.trans_matrix(1,3) = pos_x;
            this.trans_matrix(2,3) = pos_y;
            df_b_width = 1; %default blade width is 1m
            
            this.rel_pos_blade = zeros(3,df_b_width/(0.5*dx)+1);
            % notice 0.5*dx to make sure there are enough blade_pos data
            % available so that there is no skipping when diagonal
            this.rel_pos_blade(1,:) = -df_b_width/2:this.dx/2:df_b_width/2;
            % it will be -0.5 -0.45 .. ~ 0.5 for dx = 0.1
            this.rel_pos_blade(2,:) = 0.2*ones(size(this.rel_pos_blade(1,:)));
            % [x1 x2 x3 ...
            %  y1 y2 y3 ...
            %  1  1  1  ...]
            this.rel_pos_blade(3,:) = ones(size(this.rel_pos_blade(3,:)));
            
            this.rel_pos_blade_front = this.rel_pos_blade + ...
                                   [0;dx;0]*ones(1,length(this.rel_pos_blade));  
            this.rel_pos_blade_back = this.rel_pos_blade + ...
                                   [0;-dx;0]*ones(1,length(this.rel_pos_blade));                   
            
        end
        
        function move(this,forward,turn)
            % forward is 1  backward is -1 
            % turn 1 for ccw 10 deg, -1 for cw 10 deg
            this.body_ori = this.body_ori + turn*10*pi/180;
            th = turn*10*pi/180;
            p = forward*this.dx;
            
            R = [cos(th) -sin(th) 0;
                 sin(th) cos(th)  p;
                 0       0        1];
            
            this.trans_matrix = this.trans_matrix*R;
            this.body_pos(1) = this.trans_matrix(1,3);
            this.body_pos(2) = this.trans_matrix(2,3);
            
            this.blade_pos = this.trans_matrix*this.rel_pos_blade;
            this.blade_pos_front = this.trans_matrix*this.rel_pos_blade_front;
            this.blade_pos_back = this.trans_matrix*this.rel_pos_blade_back;
        end
        
        function update_pos(this, px,py,th)      
            this.body_pos = [px; py];
            
            th = th*pi/180;
            R_wb = [cos(th) -sin(th) px;
                    sin(th) cos(th)  py;
                    0       0        1] ;
 
            % blade info
            % after moving and rotating
            % find info using matrix transforms
            R_bs = zeros(3,3,length(this.blade_pos));
            for i = 1:length(this.blade_pos)
                R_bs(:,:,i) = [1 0 this.blade_pos(1,i);
                               0 1 this.blade_pos(2,i);
                               0 0 1];
            end    
            R_ws = zeros(size(R_bs));
            for i = 1:length(this.blade_pos)
                R_ws(:,:,i) = R_wb*R_bs(:,:,i);
                this.blade_pos(1,i) = R_ws(1,3,i);
                this.blade_pos(2,i) = R_ws(2,3,i);     
            end

            % blade back info
            R_bs3 = zeros(3,3,length(this.blade_pos_back));
            for i = 1:length(this.blade_pos_back)
                R_bs3(:,:,i) = [1 0 this.blade_pos_back(1,i);
                               0 1 this.blade_pos_back(2,i);
                               0 0 1];
            end    
            R_ws3 = zeros(size(R_bs3));
            for i = 1:length(this.blade_pos)
                R_ws3(:,:,i) = R_wb*R_bs3(:,:,i);
                this.blade_pos_front(1,i) = R_ws3(1,3,i);
                this.blade_pos_front(2,i) = R_ws3(2,3,i);     
            end              
            
            % blade front info
            R_bs2 = zeros(3,3,length(this.blade_pos_front));
            for i = 1:length(this.blade_pos_front)
                R_bs2(:,:,i) = [1 0 this.blade_pos_front(1,i);
                               0 1 this.blade_pos_front(2,i);
                               0 0 1];
            end    
            R_ws2 = zeros(size(R_bs2));
            for i = 1:length(this.blade_pos)
                R_ws2(:,:,i) = R_wb*R_bs2(:,:,i);
                this.blade_pos_front(1,i) = R_ws2(1,3,i);
                this.blade_pos_front(2,i) = R_ws2(2,3,i);     
            end            
            
            
        end
        
        
    end
end