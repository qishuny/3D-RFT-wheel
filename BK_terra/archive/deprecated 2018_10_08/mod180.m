function modangle = mod180(angle)

modangle = mod(angle, 360);
modangle(modangle > 180) = modangle(modangle > 180) - 360; 

end