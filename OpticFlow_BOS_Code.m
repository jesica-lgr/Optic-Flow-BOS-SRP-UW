function [IRefraction, DensityTest, Rmin_matrix] = BOS(u,v)

tic 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%BOS Algorithm%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%UW SRP Values%%%%%%%%%%%%%%%
%Ingrese el valor de la escala en pixeles / el valor de la escala en m
% esc = 5883.7;
esc = 14488.18;             %small nozzle 
%esc = 14173.228;          %big nozzle high - pressure

%esc = 12283.46457;         %unstable case with forebody
% 
% h=0.004;            		%Size of the schileren's object in meters
h = 0.0020828;

% L=0.52;                     %Distance between schlieren object and background
L = 0.24;

% Mag = 0.058837;   %Magnification Value
%Mag = 0.144881;
%Mag = 0.141173;
Mag = 0.1228346457;


%n_0=1.000297;  %Ingrese el valor del indice de refracción del medio que rodea el fenómeno 1.0002921
n_0=1.0002279;
rho_0=1.25;  %Kg/m^3		Ingrese el valor de la densidad del medio que rodea al fenomeno 
%rho_0= 7.057;  

% G = 5.3236e-4;               %Gladstone-DAle Constant of butane
G = 0.000238;                   %Nitrogen
%G = 0.0002985; 

Matxx=u;
Matyy=v; %Esto es sólo para probar que funciona Poisson con PIV

% Matxx = data_PIV_x;
% Matyy = data_PIV_y;

m=2;
a=1;
oo=1;
Dm=DifFinFun_2(Matxx,Matyy,m,a)/esc;

lx=size(Dm,1);
ly=size(Dm,2);

k= (n_0)/(G*Mag*L*h);
RHS=real((k)*Dm);

% [Density]=PoissonDN_1(real(RHS),lx,ly,m,rho_0,oo,a);
%[DensityTest]=Poisson(real(RHS),lx,ly,m,rho_0,oo,a);
[DensityTest, Rmin_matrix]=PoissonAll(real(RHS),lx,ly,m,rho_0,oo,a);

k_n= (n_0)/(Mag*L*h);
RHS_n=real((k_n)*Dm);
% [IdR]=PoissonDN_1(real(RHS_n),lx,ly,m,rho_0,oo,a);
IRefraction=DensityTest*G+1;

scale_inv = 1/esc;

step = 3;
%gráfica
figure11 = figure;
axes1 = axes('Parent',figure11);
% hold(axes1,'on');
%surf(Densidad)
mesh_metersUW(DensityTest,step,scale_inv)
shading interp
% s.EdgeColor = 'none'
view(0,-90)
%colorbar('peer',axes1);


figure22 = figure;
axes2 = axes('Parent',figure22);
% hold(axes2,'on');
%surf(IRefraccion)
mesh_metersUW(IRefraction,step,scale_inv)
shading interp
% s.EdgeColor = 'none'
view(0,-90)
%colorbar('peer',axes2);



toc




