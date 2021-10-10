%%
classdef TrajectoryContext < handle
    
    properties
        
    end
    %%
    methods
        % constructor synchronize with hmapCellContext at the beginning
        function this = TrajectoryContext()
        end
        
        function [X,Y,Z,dX,dY,dZ,theta] = minJerk(this,initial_pos,initial_vel,initial_acc,...
                                                    final_pos,final_vel,final_acc,tfinal,dT) %pos0 = [0; 0; 0];
            % pos0 vel0 acc0, posf velf accf
            xcoeff = minJerkCoefficients(initial_pos(1),initial_vel(1),initial_acc(1),final_pos(1),final_vel(1),final_acc(1),tfinal)';
            ycoeff = minJerkCoefficients(initial_pos(2),initial_vel(2),initial_acc(2),final_pos(2),final_vel(2),final_acc(2),tfinal)';
            zcoeff = minJerkCoefficients(initial_pos(3),initial_vel(3),initial_acc(3),final_pos(3),final_vel(3),final_acc(3),tfinal)';
            
%             dT = 0.05;
            X = []; % todo: allocate X,Y,Z,T 
            Y = [];
            Z = [];
            for t = 0:dT:tfinal
                T = [1; t; t^2; t^3; t^4; t^5];
                X = [ X; xcoeff*T];
                Y = [ Y; ycoeff*T];
                Z = [ Z; zcoeff*T];
            end
            % figure
            % plot3(X,Y,Z)
            dX = [];
            dY = [];
            dZ = [];
            for t = 0:dT:tfinal
                T = [0; 1; 2*t; 3*t^2; 4*t^3; 5*t^4];
                dX = [ dX; xcoeff*T];
                dY = [ dY; ycoeff*T];
                dZ = [ dZ; zcoeff*T];
            end
            
            theta = zeros(size(dX));
            theta(1) = atan2(dY(1),dX(1))-pi/2;
            for i = 2:length(dX)
                theta(i) = atan2(dY(i),dX(i))-pi/2;
                if abs(theta(i)-theta(i-1))>pi
                    theta(i) = theta(i) + 2*pi;
                end
            end
        end
        


    end
    
  
end