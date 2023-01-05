clc
clear all
close all

tic

S = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/Stable CaseCropped';
P = dir(S); %P = 'absolute or relative path to the folder S'
% mm = 266;
% nn = 171;

mm = 89;
nn = 57;

d = dir([P(1).folder '/*.tif']);
d = {d.name};
numpics = length(d);
realnum = numpics/2;
oldFolder = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode';

if numpics == 0
    disp('No suitable files found')
    cd(oldFolder)
    return
end


data_PIV_x = load('U_stable.mat'); %load u and v DaVis displacement matrixes for all frames
data_PIV_y = load('W_stable.mat');

data_PIV_x = data_PIV_x.U;
data_PIV_y = data_PIV_y.W;

%Initialize arrays
% ImgSeq =  zeros(mm,nn,2);
% MATX = zeros(mm,nn,realnum);
% MATY = zeros(mm,nn,realnum);
% U_matrix =  zeros(mm,nn,realnum);
% V_matrix =  zeros(mm,nn,realnum);
% Density_matrix=  zeros(mm-3,nn-3,realnum);
% RefractiveI_matrix =  zeros(mm-3,nn-3,realnum);

U_matrix_DaVis =  zeros(mm,nn,realnum);
V_matrix_DaVis =  zeros(mm,nn,realnum);
Density_matrix_DaVis=  zeros(mm-3,nn-3,realnum);
RefractiveI_matrix_DaVis =  zeros(mm-3,nn-3,realnum);

window = 6; %for SRP images
step = 3;

limit = 1;
diff_limit_lower = 0.05;
diff_limit_upper = 0.5;


for pict=4:(realnum)
    
%     framenum =pict-1;
%     fprintf('\tReading frame %d of %d...\n',framenum,realnum);
%     
%     left = (pict*2)-1;
%     right = pict*2;
%     
%     dleft  = d{left};
%     dright  = d{right};
%         
%     frame1 = imread(d{left});
%     frame2 = imread(d{right});
%         
%     ImgSeq(:,:,1) = double(frame1);
%     ImgSeq(:,:,2) = double(frame2);
%     
%     cd(oldFolder);
%     
%     lx=1;
%     ly=1;
% 
%     opt.eta = 0.1;
%     [Matx, Maty] = estimateOpticFlow2DUW(ImgSeq,opt);
%     
%     MATX(:,:,pict) = Matx; %pair of frames number pict
%     MATY(:,:,pict) = Matx;
%     
%     n = size(data_PIV_y(:,:,pict),2);%number of columns of each DaVis matrix
% 
%     u = zeros(size(Matx(:,:))); %Initial u, v matrixes for each frame 
%     v = zeros(size(Maty(:,:)));
%     
%     for ii = 1:n
%         %Define limits of w - columns of OF data
%         left_s = (ii*step) - 2;
%         right_s = (ii*step);
%     
%     %     left_s = (ii-1)*step + 1;
%     %     right_s = (ii*step);
%     
%         for w = left_s:right_s
%             [u(:,w),v(:,w)] = Scaling_mean_DaVisOF(Matx(:,w),Maty(:,w), data_PIV_x(:,ii,pict),data_PIV_y(:,ii,pict), limit, diff_limit_lower, diff_limit_upper,window,step);
%         end   
%     end
%     
%     U_matrix(:,:,pict) = u; %pair of frames number pict
%     V_matrix(:,:,pict) = v;

    U_matrix_DaVis(:,:,pict) = data_PIV_x(:,:,pict); %pair of frames number pict
    V_matrix_DaVis(:,:,pict) = data_PIV_y(:,:,pict);
    
    
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


    n_0=1.0002921;              %Ingrese el valor del indice de refracción del medio que rodea el fenómeno 1.0002921
    rho_0=1.204;  %Kg/m^3		Ingrese el valor de la densidad del medio que rodea al fenomeno 
    %            
    % G = 5.3236e-4;               %Gladstone-DAle Constant of butane
    G = 0.000238;                   %Nitrogen


%     Matxx=u;
%     Matyy=v; %Esto es sólo para probar que funciona Poisson con PIV

    Matxx = data_PIV_x(:,:,pict);
    Matyy = data_PIV_y(:,:,pict);

    m=2;
    a=1;
    oo=1;
    Dm=DifFinFun_2(Matxx,Matyy,m,a)/esc;

    lx=size(Dm,1);
    ly=size(Dm,2);

    k= (n_0)/(G*Mag*L*h);
    RHS=real((k)*Dm);

    % [Density]=PoissonDN_1(real(RHS),lx,ly,m,rho_0,oo,a);
%     [Density]=Poisson(real(RHS),lx,ly,m,rho_0,oo,a);
    [Density]=PoissonAll(real(RHS),lx,ly,m,rho_0,oo,a);
    Density_matrix_DaVis(:,:,pict) = Density; %pair of frames number pict

    k_n= (n_0)/(Mag*L*h);
    RHS_n=real((k_n)*Dm);
    % [IdR]=PoissonDN_1(real(RHS_n),lx,ly,m,rho_0,oo,a);
    RefractiveI_matrix_DaVis(:,:,pict)=Density*G+1;

    scale_inv = 1/esc;

%     step = 3;
%     %gráfica
%     figure11 = figure;
%     axes1 = axes('Parent',figure11);
%     % hold(axes1,'on');
%     %surf(Densidad)
%     %mesh_metersUW(Density,step,scale_inv)
%     mesh_metersUW(Density_matrix(:,:,pict),step,scale_inv)
%     shading interp
%     % s.EdgeColor = 'none'
%     view(0,-90)
%     %colorbar('peer',axes1);


%     figure22 = figure;
%     axes2 = axes('Parent',figure22);
%     % hold(axes2,'on');
%     %surf(IRefraccion)
%     %mesh_metersUW(RefractiveI,step,scale_inv)
%     mesh_metersUW(RefractiveI_matrix(:,:,pict),step,scale_inv)
%     shading interp
%     % s.EdgeColor = 'none'
%     view(0,-90)
%     %colorbar('peer',axes2);
  
end

filenameU_matrix = 'C:\Users\Jesica González\Documents\School\UNAM\Titulación\Thesis\UW Research Stay\BOSCode\U_matrix_DaVis.mat\\';
filenameV_matrix = 'C:\Users\Jesica González\Documents\School\UNAM\Titulación\Thesis\UW Research Stay\BOSCode\V_matrix_DaVis.mat\\';
filenameDensity_matrix = 'C:\Users\Jesica González\Documents\School\UNAM\Titulación\Thesis\UW Research Stay\BOSCode\Density_matrix_DaVis.mat\\';
filenameRefractiveI_matrix = 'C:\Users\Jesica González\Documents\School\UNAM\Titulación\Thesis\UW Research Stay\BOSCode\RefractiveI_matrix_DaVis.mat\\';
save(filenameU_matrix,'U_matrix_DaVis');
save(filenameV_matrix,'V_matrix_DaVis');
save(filenameU_matrix,'Density_matrix_DaVis');
save(filenameV_matrix,'RefractiveI_matrix_DaVis');

toc


