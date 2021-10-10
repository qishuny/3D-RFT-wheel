%%
function [gridX, gridY] =  makegrid(minx,maxx,miny,maxy,n)
    terrain_sizex = maxx-minx;
    xmin = minx;
    xmax = maxx;
    gridx = linspace(xmin,xmax,n*terrain_sizex+1)';
    gridX = [];
    for i = 1:terrain_sizex*n+1
        gridX = [gridX gridx];
    end

    terrain_sizey = maxy-miny;
    ymin = miny;
    ymax = maxy;
    gridy = linspace(ymin,ymax,n*terrain_sizey+1);
    gridY = [];
    for i = 1:terrain_sizey*n+1
        gridY = [gridY; gridy];
    end
end