%%
classdef SimulatorContext < handle
    % Simulator context includes pointer to everything and steps forward
    % in simulation for a given input from user or agent
    properties
        hmapPtr; 
        indexHandlerPtr;
        bodyPosePtr;
    end
    methods
        function this = SimulatorContext(hmap, indexhandler, bodypose)
            this.hmapPtr = hmap; 
            this.indexHandlerPtr = indexhandler;
            this.bodyPosePtr = bodypose;
        end
        
        function simulateStep(this, th, y)
            x0 = this.bodyPosePtr.get_x();
            y0 = this.bodyPosePtr.get_y();

            moveHT = get_HT4(0., y, 0., th);
            % postmultiply (move coord by y amount and th amount)
            this.bodyPosePtr.HT4 = this.bodyPosePtr.HT4 * moveHT;

            xf = this.bodyPosePtr.get_x();
            yf = this.bodyPosePtr.get_y();
            % points of collision in world coord
            pWorld = this.bodyPosePtr.HT4 * this.bodyPosePtr.child.HT4 * this.bodyPosePtr.child.particles;
            % non-repeating indices of collision point indices
            [indices, zvec] = this.indexHandlerPtr.getIndex(pWorld);
            dX = xf - x0;
            dY = yf - y0;
            % if it is rotating, use oreintation to move sand infront of
            % the blade
            if (y == 0)
                th = this.bodyPosePtr.get_th();
                dX = -sin(th);
                dY = cos(th);
            end
            velocity = [dX, dY];

            this.hmapPtr.resolve_collision(indices,zvec,velocity);
            this.hmapPtr.Steadystate(indices, zvec, velocity);
        end
                    
    end
end

