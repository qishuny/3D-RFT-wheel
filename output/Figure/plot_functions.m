% Define our data.
theta = linspace(0,(2*pi),30);
x=cos(theta);
y=sin(theta);
energy = linspace(0,1000,30);
numPoints = length(x);
% Get a colormap, a unique color for every energy level
cmap = jet(numPoints); % Initialize jet colormap.
% Get energy in the range 1 to numPoints so we can use that to get a row from the colormap.
qEnergy = imquantize(energy, numPoints);
for k = 1 : numPoints
	% Get the color for this energy level:
	thisEnergy = qEnergy(k);
	thisColor = cmap(thisEnergy);
	fprintf('Plotting point #%d at (%.3f, %.3f) with color (%.3f, %.3f, %.3f)\n',...
		k, x(k), y(k), cmap(k, 1), cmap(k, 2), cmap(k, 3));
	plot(x(k), y(k), '.', 'Color', cmap(k, :), 'MarkerSize', 40);
	hold on;
end
grid on;
axis square;
fprintf('Done running %s.m ...\n', mfilename);