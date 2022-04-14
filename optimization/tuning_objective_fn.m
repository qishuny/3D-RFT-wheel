function error = tuning_objective_fn(wheeldata, x, all_results)
%all_results is the data to be tuned over, if you want to tune a and b for
%a subset of data you must crop it down before passing it to this function

%perform the params setup before passing it to this function - make sure to
%set rover and soil properties, but then nullify the values for a0, a1, b0
%and b1

%initialize the error to 0
error = 0;

sf1 = x;
% sf2 = x(2);

n = length(all_results);
for i = 1:n
    %Extract measured forces to compare model prediction to
   
    result = all_results(i);
    wr = result.Vry;
    radius = 62.5;
    w = wr/radius;
    sinkage = abs(result.avg_Z);
    slipAngle = result.beta * pi / 180;
    avg_Fx = -result.avg_Fx;
    avg_Fy = -result.avg_Fy;
    avg_Fz = -result.avg_Fz;
    
    %0.005 sf1 0.04 sf2
    
    %Compute forces
%     Fx = 0;
%     Fy = 0;
%     Fz = 0;
%     [Fx, Fy, Fz] = RFT3DDEMfunc(wheeldata, slipAngle, wr, sinkage, sf1, sf2);
    
    [forces] = RFTSandTuning(wheeldata, radius, slipAngle, w, 10, sinkage, wr, 0.6, sf1, 0);
    
%     err = (Fx - avg_Fx)^2/(avg_Fx^2) + (Fy - avg_Fy)^2/(avg_Fy^2) + (Fz - avg_Fz)^2/(avg_Fz^2);
    if result.beta <= 60
        err = (forces(1) - avg_Fx)^2 + (forces(2) - avg_Fy)^2;
        error = error + err;
    end
    
    
end
error = double(error);    
end

%x - [.6 .2 -.2 -.3 0 0]';


    
    