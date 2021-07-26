clear
close all


% try bunch of stuff here
% make Heightmap, Bulldozer,
% change stuff move stuff if my classes are working properly


%% visualization practice
figure(1)
resizeFig(1.5,1.5)

% surfs = {};
% setseed(123);
% surfs{end+1} = surface; 


% % terrain = surface;
% terrain_fh = @terrain;  % working example based on wmrde
% surfs = feval(terrain_fh);
% drawSurfaces(surfs); %,anim.h_axis);

% hmap = zeros(10*4+1);
% % terrain_f = @BK_terrain;
% % surfs = feval(terrain_f);  %why @BK_terrain and feval again??
% surfs = BK_terrain(hmap);
% drawSurfaces(surfs);
% 
% makeLegible(14)
% axis off
% 
% %     view(0,0)
% view(30,30)
% zoom(1.5)


Hmap = Heightmap();
Hmap.set_param(30,0,4);
Hmap.matrix(1,1) = -0.5;
% Bull = Bulldozer();
% Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front,Bull.blade_pos_back);

load('CMUdata.mat')
for i = 1:length(CMUdata)
    Hmap.add_sand(CMUdata(i,1),CMUdata(i,2),2);
end
Hmap.visualize;
colormap Copper;

view(30,30)
zoom(1.5)
Hmap.visualize_dynamics();

% hmap = Hmap.matrix;
% surfs = BK_terrain(hmap);
% figure()
% drawSurfaces(surfs);


% 
% terrain_sizex = 4; %m
% terrain_sizey = 4; %m
% 
% % % Height Map Initialization
% n = 10; % grid per m
% dx = 1/n; % grid size note to self: 1/25 fails 
% HeightMap = zeros(terrain_sizex*n+1,terrain_sizey*n+1);
% HeightMap_comp = zeros(size(HeightMap));
% [wid, leng] = size(HeightMap);
% 
% 
% % % Some Initial Sand Piled up
% ini_vol = 0.8;
% ini_h = ini_vol/(dx^2);
% for i = 20:30
%     HeightMap(i,20) = ini_h/8;
% end
% for j = 20:4
%     HeightMap(20,j) = ini_h/8;
% end
% HeightMap;
% 
% % for i = 1:2
% %     for j = 1:2
% %     HeightMap(int8(wid/2)+(-1)^i,int8(leng/2)+(-1)^j) = ini_h/5;
% %     end
% % end
% % max(HeightMap)
% 
% % Setting up the Grid (x,y location) for my Height Map Visualization
% gridx = linspace(0,terrain_sizex,n*terrain_sizex+1)';
% gridX = [];
% for i = 1:terrain_sizex*n+1
%     gridX = [gridX gridx];
% end
% gridy = linspace(0,terrain_sizex,n*terrain_sizey+1);
% gridY = [];
% for i = 1:terrain_sizey*n+1
%     gridY = [gridY; gridy];
% end

%% % ----------------- Discrete Element Simulation ------------------------
% % Just One big Block
% my_block = DiscreteElement();
% my_block.assign_attributes(0.1,1000,[1,1,1]);
% 
% my_block_vel = [1 0 0]'; %say Im pushing with buldozer at constant speed
% 
% 
% for i = 1:10
%     my_block.update(my_block_vel,0.1);
%     my_block.position(1)
%     pos = my_block.position(1:2);
%     h = my_block.position(3);
%     
%     HeightMap = zeros(terrain_sizex*n+1,terrain_sizey*n+1);
%     HeightMap(int8(pos(1)/dx),int8(pos(2)/dx)) = h;
%     surf(gridx,gridy,HeightMap)
%     axis([0 terrain_sizex 0 terrain_sizey 0 terrain_sizex])
%     drawnow
% end





%% % ----------------- Compression Simulation ---------------
% for i = 1:10
%     position = {[20+i 20] [20+i 21] [20+i 22] [20+i 23]};
%     [update updatecomp] = update_height_compression(HeightMap,HeightMap_comp, 1000, position, dx);
%     HeightMap = HeightMap + update;
%     HeightMap_comp = HeightMap_comp + update;
%     surf(gridx,gridy,HeightMap)
%     axis([0 terrain_sizex 0 terrain_sizey -0.1 0.1])
%     drawnow
% end
% 
% for i = 1:10
%     position = {[30-i 17] [30-i 18] [30-i 19] [30-i 20]};
%     [update updatecomp] = update_height_compression(HeightMap,HeightMap_comp, 1000, position, dx);
%     HeightMap = HeightMap + update;
%     HeightMap_comp = HeightMap_comp + update;
%     surf(gridx,gridy,HeightMap)
%     axis([0 terrain_sizex 0 terrain_sizey -0.1 0.1])
%     drawnow
% end


%% % ----------------- Sand Slipping Simulation -------------------
% each for loop is essentially time increment
% for i = 1:150
%     
%     update_Height = update_height(HeightMap,dx);
%     HeightMap = HeightMap + update_Height;
%     surf(gridx,gridy,HeightMap)
%     axis([0 terrain_sizex 0 terrain_sizey 0 terrain_sizex])
%     drawnow
% end
% 
% gridMap = zeros(size(HeightMap));
% 
% surf(gridx,gridy,HeightMap)
% axis ([0 terrain_sizex 0 terrain_sizey 0 terrain_sizex])
% max(max(HeightMap))


%% ----- Sand Simulation with Smaller Grid ------

% Hmap = Heightmap();
% Hmap.set_param(40,0,4);
% size(Hmap.matrix)
% Bull = Bulldozer();
% Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front);
% 
% for i = 1:10
%     Hmap.add_sand(1.4,1.4+i*0.1,4);
% end
% Hmap.visualize();
% Hmap.matrix;
% pause
% for i = 1:10
% Hmap.visualize_onestep();
% pause
% end

% Hmap.visualize_dynamics();


% load('CMUdata.mat')
% for i = 1:length(CMUdata)
%     Hmap.add_sand(CMUdata(i,1),CMUdata(i,2),2);
% end
% Hmap.visualize_dynamics();


%% ---- Checking Near Steady State ----
Hmap = Heightmap();
volume = 0.1;
n = 50;
dx = 1/n;
height = volume/dx^2;
Hmap.set_param(n,0,1.5);
Hmap.add_sand(0.75,0.75,height);
% Hmap.steady_state;
% Hmap.visualize
for i = 1:2000
Hmap.visualize_onestep(0.01);
end
% Hmap.visualize_dynamics()
% max(max(Hmap.matrix))


% for i = 1:200
%     Hmap.visualize_onestep()
%     pause
% end
% Hmap.visualize_dynamics();





%% ----- Sand Simulation CMU!!! -----
Hmap = Heightmap();
Hmap.set_param(30,0,4);
% Bull = Bulldozer();
% Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front,Bull.blade_pos_back);

% for i = 1:10
%     Hmap.add_sand(1.4,1.4+i*0.1,20);
% end
% Hmap.visualize_dynamics();

load('CMUdata.mat')
for i = 1:length(CMUdata)
    Hmap.add_sand(CMUdata(i,1),CMUdata(i,2),2);
end
Hmap.visualize;
% pause
% Hmap.visualize_onestep();
% pause
Hmap.visualize_dynamics();


%% ------------- Impulse Response --------------------
clear
close all

Hmap = Heightmap();
Hmap.set_param(40,0,2);
volume = 0.08;
impulse_w = 0.2;
impulse_l = 0.2;
impulse_h = volume/(impulse_w*impulse_l)

x_idx = (impulse_w / Hmap.dx);
y_idx = (impulse_l / Hmap.dx);

for i = 1:x_idx
    for j = 1:y_idx
        Hmap.add_sand(0.9+Hmap.dx*(i-1),0.9+Hmap.dx*(j-1),impulse_h);
    end
end

% Hmap.add_sand(2,2,impulse_h);
Hmap.visualize
axis([0.5 1.5 0.5 1.5 -0.1 1.4])
% view([0 0])
pause

tic
Hmap.update(10)
toc
Hmap.visualize
axis([0.5 1.5 0.5 1.5 -0.1 2])
axis([0.5 1.5 0.5 1.5 -0.1 1.4])
% view([0 0])
pause

tic
Hmap.update(20)
toc
Hmap.visualize
axis([0.5 1.5 0.5 1.5 -0.1 2])
axis([0.5 1.5 0.5 1.5 -0.1 1.4])
% view([0 0])
pause

% Hmap.visualize_dynamics
tic
Hmap.steady_state
toc
Hmap.visualize
axis([0.5 1.5 0.5 1.5 -0.1 2])
axis([0.5 1.5 0.5 1.5 -0.1 1.4])
max(max(Hmap.matrix))
% view([0 0])
pause 

tic
Hmap.steady_state
Hmap.steady_state
toc
Hmap.visualize
axis([0.5 1.5 0.5 1.5 -0.1 2])
axis([0.5 1.5 0.5 1.5 -0.1 1.4])
max(max(Hmap.matrix))
% view([0 0])

for i = 1:10
Hmap.steady_state
end
Hmap.visualize
axis([0.5 1.5 0.5 1.5 -0.1 2])
axis([0.5 1.5 0.5 1.5 -0.1 1.4])
max(max(Hmap.matrix))
% view([0 0])
%% -------------- sand pushing with parts -----------------
clear
close all

Hmap = Heightmap();
Hmap.set_param(40,0,4);
Bull = Bulldozer();
blade_z = -0.1;

Bull.parts{1}.name = "wheel1";
Bull.parts{1}.pos = [-0.2121 0.1414]; %x, y pos rel to vehicle
Bull.parts{1}.length = 0.8;
Bull.parts{1}.orientation = 45; % in degrees
Bull.parts{1}.z = -0.1;

Bull.parts{2}.name = "wheel2";
Bull.parts{2}.pos = [0.2121 0.1414]; %x, y pos rel to vehicle
Bull.parts{2}.length = 0.8;
Bull.parts{2}.orientation = -45; % in degrees
Bull.parts{2}.z = -0.1;

Bull.synchronize(Hmap.dx,1.5,0.5);

for i = 1:20
    Bull.move(1,0);
    Bull.pos
%     Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front,Bull.blade_pos_back,blade_z);
    Hmap.synchronize(Bull.rw_blade,Bull.rw_blade_f,Bull.rw_blade_b,blade_z);
%     Bull.rw_blade
%     Bull.rw_
    Hmap.steady_state();
    Hmap.visualize();
% %     Hmap.visualize_dynamics
% pause
end

% for i = 1:3
%     Bull.move(0,-1);
%     Bull.pos
% %     Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front,Bull.blade_pos_back,blade_z);
%     Hmap.synchronize(Bull.rw_blade,Bull.rw_blade_f,Bull.rw_blade_b,blade_z);
%     Hmap.steady_state();
%     Hmap.visualize();
% % %     Hmap.visualize_dynamics
% end
% 
% for i = 1:20
%     Bull.move(1,0);
%     Bull.pos
% %     Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front,Bull.blade_pos_back,blade_z);
%     Hmap.synchronize(Bull.rw_blade,Bull.rw_blade_f,Bull.rw_blade_b,blade_z);
%     Hmap.steady_state();
%     Hmap.visualize();
% % %     Hmap.visualize_dynamics
% end


%% -------------- Sand Simulation Bulldozer Pushing ----------------------------
clear
close all

Hmap = Heightmap();
Hmap.set_param(40,0,4);
Bull = Bulldozer();
blade_z = -0.1;
Bull.synchronize(Hmap.dx,1.5,0.5);
% Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front,Bull.blade_pos_back,blade_d);

% for i = 1:9
%     Bull.move(0,1);
%     Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front,Bull.blade_pos_back,blade_d);
%     
% end
% Hmap.visualize
% Bull.body_ori
% Bull.body_pos

volume = 0.2;
impulse_w = 0.2;
impulse_l = 0.2;
impulse_h = volume/(impulse_w*impulse_l);

x_idx = (impulse_w / Hmap.dx);
y_idx = (impulse_l / Hmap.dx);

for i = 1:x_idx
    for j = 1:y_idx
        Hmap.add_sand(1.9+Hmap.dx*(i-1),1.9+Hmap.dx*(j-1),impulse_h);
    end
end

for i = 1:x_idx
    for j = 1:y_idx
        Hmap.add_sand(1.2+Hmap.dx*(i-1),3+Hmap.dx*(j-1),impulse_h);
    end
end

for i = 1:x_idx
    for j = 1:y_idx
        Hmap.add_sand(2.6+Hmap.dx*(i-1),1.5+Hmap.dx*(j-1),impulse_h);
    end
end
for i = 1:40
        Hmap.add_sand(0.2+Hmap.dx*(i-1),2+Hmap.dx*(j-1),-impulse_h/5);
end

Hmap.visualize_dynamics();
% tic
% Hmap.steady_state();
% toc
% Hmap.visualize


% for i = 1:10
%     Hmap.add_sand(3,1+i*0.1,20);
% end
% Hmap.visualize_dynamics();
% Hmap.visualize_dynamics();

% tic
for i = 1:10
    Bull.move(1,0);
    Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front,Bull.blade_pos_back,blade_z);
    Hmap.steady_state();
    Hmap.visualize();
% %     Hmap.visualize_dynamics
end
% pause

for i = 1:3
    Bull.move(0,-1);
    Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front,Bull.blade_pos_back,blade_z);
    Hmap.steady_state
    Hmap.visualize
%     pause
end

for i = 1:35
    Bull.move(1,0);
    Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front,Bull.blade_pos_back,blade_z);
    Hmap.steady_state();
    Hmap.visualize();
end
% pause
for i = 1:35
    Bull.move(1,0);
    Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front,Bull.blade_pos_back,blade_z);
    Hmap.steady_state();
    Hmap.visualize();
end
% pause

for i = 1:18
    Bull.move(0,1);
    Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front,Bull.blade_pos_back,blade_z);
    Hmap.steady_state
    Hmap.visualize
%     pause
end


for i = 1:70
    Bull.move(1,0);
    Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front,Bull.blade_pos_back,blade_z);
    Hmap.steady_state();
    Hmap.visualize();
end
% toc
% Hmap.steady_state
% Hmap.visualize


% 
% % for i = 1:10
% %     Hmap.add_sand(1.2,1.4+i*0.1,20);
% % end
% % Hmap.visualize_dynamics();
% 
% % lets make ridge
% % for i = 1:10
% %     Hmap.add_sand(1.4,1.4+i*0.1,20);
% % end
% % Hmap.visualize_dynamics();
% 
% % for i = 1:10
% %     Bull.move(1,0);
% %     Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front);
% %     Hmap.visualize_dynamics();
% % end
% 
% % Bull.move(0,-3);
% 
% 
% 

% 
% for i = 1:10
%     Bull.move(1,0);
%     Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front);
%     Hmap.visualize_dynamics();
% end
% Bull.move(1,0);
% Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front);
% Hmap.visualize_dynamics();
% Bull.move(1,0);
% Hmap.synchronize(Bull.blade_pos,Bull.blade_pos_front);
% Hmap.visualize_dynamics();

% for i = 1:10
%     Hmap.add_sand(2,2.4+i*0.1,1);
% end
% Hmap.visualize_dynamics();
% 
% for i = 1:10
%     Hmap.add_sand(2,2.6+i*0.1,1);
% end
% Hmap.visualize_dynamics();
% 
% for i = 1:10
%     Hmap.add_sand(2,2.6+i*0.1,1);
% end
% Hmap.visualize_dynamics();
% 
% for i = 1:10
%     Hmap.add_sand(2,2.6+i*0.1,3);
% end
% Hmap.visualize_dynamics();


%% % --------------- Sand Slipping with Matrix Method ---------------------
% blade_idx = [18 18 18 18 18 18;20 21 22 23 24 25];
% for i = 1:150
% %     HeightMap = mtx_method_update_height(HeightMap,dx);
%     HeightMap = mtx_method_update_height(HeightMap,dx,blade_idx);
%     surf(gridX,gridY,HeightMap)
%     axis([0 terrain_sizex 0 terrain_sizey 0 terrain_sizex])
%     drawnow  
% 
% end

% prompt = 'What is the original value? ';
% x = input(prompt)
% y = x*10


% %% %----------------- Bulldozer Simulation ----------------------------
% % % blade push info, blade height info
% 
% data = load('HeightMap_4040.mat');
% HeightMap = data.HeightMap;
% 
% blade_info = cell(1,10);
% for i = 1:length(blade_info)
% blade_info{i} = [11 10+i-1; 12 11+i-1];
% end
% blade_h = 0;
% 
% for bulldozer_move = 1:20
%     % bulldozer moved
%     for i = 1:length(blade_info)
%     blade_info{i} = [13+bulldozer_move 10+i-1; 14+bulldozer_move 11+i-1];
%     end
%     blade_h = -0.1;
%     
%     % bulldozer update
%     for i = 1:length(blade_info)
%     dH = bulldozer_update(HeightMap,blade_info, blade_h);
%     end
%     HeightMap = HeightMap + dH;
%     surf(gridx,gridy,HeightMap)
%     axis ([0 terrain_sizex 0 terrain_sizey -0.2 terrain_sizex])
% 
%     % heightMap update
%     max_hdiff = 1;
%     while max_hdiff > 0.01
%         update_Height = update_height(HeightMap,dx,blade_info);
%         HeightMap = HeightMap + update_Height;
%         max_hdiff = max(max(abs(update_Height)))
%     end
%     
%     surf(gridx,gridy,HeightMap) 
%     axis([0 terrain_sizex 0 terrain_sizey -0.2 terrain_sizex])
%     drawnow
%     
% %     for i = 1:50
% %         update_Height = update_height(HeightMap,dx,blade_info);
% %         HeightMap = HeightMap + update_Height;
% %         surf(gridx,gridy,HeightMap)
% %         axis([0 terrain_sizex 0 terrain_sizey -0.2 terrain_sizex])
% %         drawnow
% %     end
%     
%     
% end
