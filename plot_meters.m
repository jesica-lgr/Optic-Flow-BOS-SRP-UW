function plot_meters(data, step)

[rows, cols] = size(data);
divisions = floor(rows/step); %number of intervals/divisions in PIVLab plot

xStart = 1;
dx = 1;
N = rows;
x = xStart + (0:N-1)*dx; %list of divisions in px units

X = x*step*0.000045;%convert px to m

plot(X, data);
