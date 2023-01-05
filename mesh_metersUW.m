function mesh_metersUW(data, step,scale_inv)

[rows, cols] = size(data);
divisions = floor(rows/step); %number of intervals/divisions in PIVLab plot

xStart = 1;
dx = 1;
%N = rows;
N = cols;
x = xStart + (0:N-1)*dx; %list of divisions in px units
yStart = 1;
dy = 1;
% M = cols;
M = rows;
y = yStart + (0:M-1)*dy; %list of divisions in px units

% X = x*step*scale_inv;%convert px to m
% Y = y*step*scale_inv;%convert px to m
X = x*scale_inv;%convert px to m
Y = y*scale_inv;%convert px to m

mesh(X,Y, data);