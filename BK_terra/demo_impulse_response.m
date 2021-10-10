%% ---- Hmap Sand Rolling Demo V2 --
clear 
close all
n = 10; %grid per m
minx = -8;
miny = -8;
maxx = 8;
maxy = 8;
Sand = Hmap(minx,miny,maxx,maxy,n);
% Sand.dT = 2;

%-----------initialize visualizeContext
visualizeContext = VisualizeContext(Sand.minx,Sand.miny,Sand.maxx,Sand.maxy,Sand.n);

% ----------initialize indexHandler
indexHandler = IndexHandler(Sand.get_dx,Sand.minx,Sand.miny,Sand.maxx,Sand.maxy);

%-----------visualize hmap
% visualizeContext.visualize(Sand.matrix);

tic
% ----------perturb and visualize
% Sand.add_sand(-8,8,-20); %perturb
for i = 1:40
Sand.add_sand(i*Sand.get_dx,0,10)
end
visualizeContext.visualize(Sand.matrix);
visualizeContext.hAxes
pause
% Sand.update_onestep;
% while ~Sand.reachedSteadystate
%     Sand.update_onestep;
%     visualizeContext.visualize(Sand.matrix);
% end


%perturb
for j = 1:30
    Sand.add_sand(4,j*Sand.get_dx,10);
end
pause
visualizeContext.hAxes
visualizeContext.visualize(Sand.matrix);

%update and visualize
% Sand.update_onestep;
% while ~Sand.reachedSteadystate
%     Sand.update_onestep;
%     visualizeContext.visualize(Sand.matrix);
% end


%perturb 
for i = 1:40
    Sand.add_sand(-4+i*Sand.get_dx,-4,-5);
end
visualizeContext.visualize(Sand.matrix);
pause

%update and visualize
Sand.update_onestep;
while ~Sand.reachedSteadystate
    Sand.update_onestep;
    visualizeContext.visualize(Sand.matrix);
%     pause
end
toc