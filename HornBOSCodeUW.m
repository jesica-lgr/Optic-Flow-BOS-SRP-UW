clc
clear all
close all
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%Validation images SF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/ValidationImagesSF/A%03d.tif';


%%%%%%%%%%%%%%%%%%%%%%%Gas Of Lighter Images%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%CLOSE UP 1%%%%
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/Código Flujo óptico-BOS/Horn-BOS merge/UWLighter/105mm_CloseUp_AdvancedMovement_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/Código Flujo óptico-BOS/Horn-BOS merge/UWLighter/105mm_CloseUp_InitialMovement_%03d.tif';

%%%%CLOSE UP 2%%%%
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/Código Flujo óptico-BOS/Horn-BOS merge/UWLighter/105mm_CloseUp2_AdvancedMovement_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/Código Flujo óptico-BOS/Horn-BOS merge/UWLighter/105mm_CloseUp2_InitialMovement_%03d.tif';

%%%%Close Up 1 and great spacing between frames%%%%
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/Código Flujo óptico-BOS/Horn-BOS merge/UWLighter/105mm_CloseUp_Reference_%03d.tif';



%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UWLighter/105mm_CloseUp_ReferenceAdv_%03d.tif';


%%%%%%%%300mm Photos Lighter%%%%%%%%%%%%%%
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UWLighter/300mm_1280x1024_5micros_66cm_Lighter_Initial_%03d.tif';

% filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UWLighter/300mm_512x512_3micros_250cm_Lighter2_%03d.tif';
% filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UWLighter/300mm_512x512_3micros_250cm_Lighter1_%03d.tif';

%A bit of problems with displacements of 5px
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UWLighter/300mm_CloseUp_N2_512x512_3micros_203cm_0.7D_MaxAperture_LighterComplete_%03d.tif';

%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UWLighter/300mm_CloseUp_N2_512x512_3micros_250cm_0.5D_MaxAperture_LighterComplete_Ini_crop_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UWLighter/300mm_CloseUp_N2_512x512_3micros_250cm_0.5D_MaxAperture_LighterComplete_Adv_crop_%03d.tif';



%%%%%%%%%%%%300mm Photos SRP%%%%%%%%%%%

%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SRP/BOS_18_1.75_22.0_SF_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SRP/BOS_18_1.75_22.0_SRP_%03d.tif';


%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SRP/BOS_18_1.75_22.0_SRP_AB_%03d.tif';





%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SRP/B_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SRP/C2_%03d.tif';


%%%%%%%%%%%Changing Focus Middle, Nozzle, Pattern%%%%%%%%%
%%%%%%%%%%%%24cm and 5.1cm AB%%%%%%%%%%%%%%%%%%%%%%
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF/FM_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF/FocusMiddle_10_%03d.tif';

%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF/FocusNozzle_10_%03d.tif';

%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF/FocusPattern_10_%03d.tif';


%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF/FocusBoth_AB_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF/FocusBoth_%03d.tif';

%%%%%%%%%%%%24cm background and SF%%%%%%%%%%%%%%%%%%%%%%
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF/Focus_Middle_2_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF/Focus_Nozzle_2_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF/Focus_Pattern_2_%03d.tif';


%%%%%%%Tests 18.11.2021 SF and SRP focus Pattern%%%%%
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF/BOS_18_0.5_22.0_V2_SF_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF/BOS_18_0.5_22.0_V2_SRP_%03d.tif';
% 
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF/BOS_18_1.75_22.0_V2_SF_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF/BOS_18_1.75_22.0_V2_SRP_%03d.tif';
% 
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF/BOS_38_0.5_22.0_SF_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF/BOS_38_0.5_22.0_SRP_%03d.tif';

%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SRP/test_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SRP/08.12.2021/SRP_%03d.tif';

%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF Lens2/BOS_38_0.5_21.5_R3_%03d.tif';



%%%%%%%Tests 07.12.2021 New setup with lens 2%%%%%
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF Lens2/80mm_lens2_512x512_maxfocus_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF Lens2/135mm_lens2_800x600_maxfocus_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF Lens2/135mm_lens2_800x600_maxfocus_34.3_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SF Lens2/135mm_lens2_800x600_maxfocus_14.3_%03d.tif';

%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SRP/09.12.2021/BOS_38_0.5_21.5_R3_1_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/UW SRP/09.12.2021/BOS_38_0.5_21.5_R3_2_%03d.tif';

filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/Stable CaseCropped/Stable_Case_150_%03d.tif';




%%%%%%%%%%%%%%%%%%%%%%%Read Image Sequence%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ImgSeq = readImgSeq(filePattern,0,1);


%%%%%%%%%%%%Estimate optic flow for the image sequence%%%%%%%%%%%%%%%%%%%%%

lx=1;
ly=1;

opt.eta = 0.1;
[Matx, Maty] = estimateOpticFlow2DUW(ImgSeq,opt);

%%%%%%%%%%Scale displacement's data using PIV results%%%%%%%%%%%%%%%%%%%%%%


% data_PIV_x = load('u_PIV_validation.mat');
% data_PIV_y = load('v_PIV_validation.mat');

% data_PIV_x = load('u_Focus_Both_AB.mat');
% data_PIV_y = load('v_Focus_Both_AB.mat');

% data_PIV_x = load('u_Focus_Middle_10.mat');
% data_PIV_y = load('v_Focus_Middle_10.mat');

% data_PIV_x = load('u_Focus_Pattern_10.mat');
% data_PIV_y = load('v_Focus_Pattern_10.mat');

% data_PIV_x = load('u_Focus_Nozzle_10.mat');
% data_PIV_y = load('v_Focus_Nozzle_10.mat');


% data_PIV_x = load('u_Focus_Middle_2.mat');
% data_PIV_y = load('v_Focus_Middle_2.mat');
% data_PIV_x = load('u_Focus_Nozzle_2.mat');
% data_PIV_y = load('v_Focus_Nozzle_2.mat');
% data_PIV_x = load('u_Focus_Pattern_2.mat');
% data_PIV_y = load('v_Focus_Pattern_2.mat');

% data_PIV_x = load('u_BOS_18_0.5_22.0_V2_SF.mat');
% data_PIV_y = load('v_BOS_18_0.5_22.0_V2_SF.mat');
% data_PIV_x = load('u_BOS_18_1.75_22.0_V2_SF.mat');
% data_PIV_y = load('v_BOS_18_1.75_22.0_V2_SF.mat');
% data_PIV_x = load('u_BOS_38_0.5_22.0_SF.mat');
% data_PIV_y = load('v_BOS_38_0.5_22.0_SF.mat');

% data_PIV_x = load('u_BOS_18_0.5_22.0_V2_SRP.mat');
% data_PIV_y = load('v_BOS_18_0.5_22.0_V2_SRP.mat');
% data_PIV_x = load('u_BOS_18_1.75_22.0_V2_SRP.mat');
% data_PIV_y = load('v_BOS_18_1.75_22.0_V2_SRP.mat');
% data_PIV_x = load('u_BOS_38_0.5_22.0_SRP.mat');
% data_PIV_y = load('v_BOS_38_0.5_22.0_SRP.mat');


% data_PIV_x = load('u_BOS_38_0.5_21.5_R3_1.mat');
% data_PIV_y = load('v_BOS_38_0.5_21.5_R3_1.mat');
data_PIV_x = load('u_BOS_38_0.5_21.5_R3_2.mat');
data_PIV_y = load('v_BOS_38_0.5_21.5_R3_2.mat');




% data_PIV_x = load('u_PIV_105mm_CloseUp_InitialMovement.mat');
% data_PIV_y = load('v_PIV_105mm_CloseUp_InitialMovement.mat');

% data_PIV_x = load('u_PIV_105mm_CloseUp_AdvancedMovement.mat');
% data_PIV_y = load('v_PIV_105mm_CloseUp_AdvancedMovement.mat');
% % 
% data_PIV_x = load('u_PIV_105mm_CloseUp2_InitialMovement.mat');
% data_PIV_y = load('v_PIV_105mm_CloseUp2_InitialMovement.mat');
% 
% data_PIV_x = load('u_PIV_105mm_CloseUp2_AdvancedMovement.mat');
% data_PIV_y = load('v_PIV_105mm_CloseUp2_AdvancedMovement.mat');


% data_PIV_x = load('u_PIV_105mm_CloseUp_Reference.mat');
% data_PIV_y = load('v_PIV_105mm_CloseUp_Reference.mat');


% data_PIV_x = load('u_PIV_105mm_CloseUp_ReferenceAdv.mat');
% data_PIV_y = load('v_PIV_105mm_CloseUp_ReferenceAdv.mat'); 



data_PIV_x = data_PIV_x.u_filtered; %for Focus Middle_10 and Nozzle_10 Middle_2 and Pattern _2 BOS_18_0.5_22.0_V2_SF
data_PIV_y = data_PIV_y.v_filtered; %Also for v_BOS_38_0.5_21.5_R3_1.mat





% data_PIV_x = data_PIV_x.u_original; %For rest of images
% data_PIV_y = data_PIV_y.v_original;

% data_PIV_x = data_PIV_x.data_PIV_x;
% data_PIV_y = data_PIV_y.data_PIV_y;
data_PIV_x = data_PIV_x{1,1};
data_PIV_y = data_PIV_y{1,1};

window = 16; %for SF validation images and Nozzle_2 and Pattern_2 and all SRP
step = window/2;

% window = 16; %For gas of lighter
% step = window/2;

% window = 32; %For soplete 8bits and soplete enc and focus pattern
% step = window/2;

% window = 32; %For SF Focus Both AB and Middle_2
% step = 8;
    
% u = zeros(size(Matx,1)-10,1);
% v = zeros(size(Maty,1)-10,1);
% u = zeros(size(Matx,1)-8,1);
% v = zeros(size(Maty,1)-8,1);
%u = zeros(size(Matx,1)-(step),1);
%v = zeros(size(Maty,1)-(step),1);
% u = zeros(size(Matx,1)-(step/2),1);
% v = zeros(size(Maty,1)-(step/2),1);

n = size(data_PIV_y,2);

u = zeros(size(Matx));
v = zeros(size(Maty));

limit = 1;
diff_limit_lower = 0.05;
diff_limit_upper = 0.5;

for ii = 1:n
    %Define limits of w - columns of Horn data
    left_s = (ii*step) - ((step/2)-1);
    right_s = (ii*step) + (step/2);
    
    for w = left_s:right_s
        [u(:,w),v(:,w)] = Scaling_mean_test(Matx(:,w),Maty(:,w), data_PIV_x(:,ii),data_PIV_y(:,ii), limit, diff_limit_lower, diff_limit_upper,window,step);
    end
    
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%BOS Algorithm%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% esc=6688.5 ;                %Ingrese el valor de la escala en pixeles / el valor de la escala en m
% % 
% h=0.004;            		%Tamaño del objeto Schlieren  en metros
% L=0.0876;                	%Distancia entre objeto schlieren y el Background en metros
% % 
% Mag=0.13;               	%Ingrese el valor de la magnificación
% % 
% 
% % 
% n_0=1.0002921;              %Ingrese el valor del indice de refracción del medio que rodea el fenómeno 1.0002921
% rho_0=1.204;  %Kg/m^3		Ingrese el valor de la densidad del medio que rodea al fenomeno 
% G = 2.2644e-4;              %Consante de Gladstone-DAle

%%%%%%%%%%%%UW SRP Values%%%%%%%%%%%%%%%
%Ingrese el valor de la escala en pixeles / el valor de la escala en m
% esc = 5883.7;
%esc = 14488.18;             %small nozzle 
%esc = 14173.228;          %big nozzle high - pressure

esc = 12283.46457;         %unstable case with forebody
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


%%%%%%%%%%%%UW SF Values%%%%%%%%%%%%%%%
%Ingrese el valor de la escala en pixeles / el valor de la escala en m
%esc = 14400;
% esc = 13900;
%esc = 15100; 
%esc = 13700;


% h=0.002;            		%Size of the schileren's object in meters
% %L = 0.051;
% L=0.24;                     %Distance between schlieren object and background
% % 
% %Mag = 0.144;   %Magnification Value
% % Mag = 0.139;
% %Mag = 0.151;
% Mag = 0.137;
% 
% n_0=1.0002921;              %Ingrese el valor del indice de refracción del medio que rodea el fenómeno 1.0002921
% rho_0=1.204;  %Kg/m^3		Ingrese el valor de la densidad del medio que rodea al fenomeno 
% %            
% G = 0.000238;               %Gladstone-DAle Constant of NITROGEN





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

[Densidad]=PoissonDN_1(real(RHS),lx,ly,m,rho_0,oo,a);

k_n= (n_0)/(Mag*L*h);
RHS_n=real((k_n)*Dm);
% [IdR]=PoissonDN_1(real(RHS_n),lx,ly,m,rho_0,oo,a);
IRefraccion=Densidad*G+1;

scale_inv = 1/esc;


%gráfica
figure11 = figure;
axes1 = axes('Parent',figure11);
% hold(axes1,'on');
%surf(Densidad)
mesh_metersUW(Densidad,step,scale_inv)
shading interp
% s.EdgeColor = 'none'
view(0,-90)
%colorbar('peer',axes1);


figure22 = figure;
axes2 = axes('Parent',figure22);
% hold(axes2,'on');
%surf(IRefraccion)
mesh_metersUW(IRefraccion,step,scale_inv)
shading interp
% s.EdgeColor = 'none'
view(0,-90)
%colorbar('peer',axes2);



toc




