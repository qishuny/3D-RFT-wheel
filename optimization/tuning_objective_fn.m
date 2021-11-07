function error = tuning_objective_fn(x, all_results)
%all_results is the data to be tuned over, if you want to tune a and b for
%a subset of data you must crop it down before passing it to this function

%perform the params setup before passing it to this function - make sure to
%set rover and soil properties, but then nullify the values for a0, a1, b0
%and b1

%initialize the error to 0
error = 0;

sf1 = x(1);
sf2 = x(2);

n = length(all_results);
for i=1:n
    %Extract measured forces to compare model prediction to
   
    result = all_results(i);
    wr = result.Vry;
    sinkage = abs(result.avg_Z);
    slipAngle = result.beta * pi / 180;
    avg_Fx = -result.avg_Fy;
    avg_Fy = -result.avg_Fx;
    avg_Fz = -result.avg_Fz;
    
    %0.005 sf1 0.04 sf2
    
    %Compute forces
%     Fx = 0;
%     Fy = 0;
%     Fz = 0;
    [Fx, Fy, Fz] = RFT3DDEMfunc(slipAngle, wr, sinkage, sf1, sf2);
    

    err = (Fx - avg_Fx)^2/(avg_Fx^2) + (Fy - avg_Fy)^2/(avg_Fy^2) + (Fz - avg_Fz)^2/(avg_Fz^2);
%     err = (Fx - avg_Fx)^2 + (Fy - avg_Fy)^2 + (Fz - avg_Fz)^2;
    error = error + err;
end
error = double(error)    
end

%x - [.6 .2 -.2 -.3 0 0]';


    
    