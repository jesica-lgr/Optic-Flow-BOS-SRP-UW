clc
clear all
close all
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%Validation images SF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/ValidationImagesSF/A%03d.tif';

%%%%%%%%%%%%%%%%%%%%%%%%SRP IMAGES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%STABLE CASE
filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/UW Research Stay/BOSCode/Stable CaseCropped/Stable_Case_100_%03d.tif';


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

% data_PIV_x = data_PIV_x.u_original; %For rest of images
% data_PIV_y = data_PIV_y.v_original;

% data_PIV_x = data_PIV_x.data_PIV_x;
% data_PIV_y = data_PIV_y.data_PIV_y;
% data_PIV_x = data_PIV_x{1,1};
% data_PIV_y = data_PIV_y{1,1};

% window = 16; %for SF validation images 
% step = window/2;

data_PIV_x = load('U_stable.mat'); %load u and v DaVis displacement matrixes for all frames
data_PIV_y = load('W_stable.mat');
% 
data_PIV_x = data_PIV_x.U;
data_PIV_y = data_PIV_y.W;

V_DaVis_100 = data_PIV_y(:,:,100);
U_DaVis_100 = data_PIV_x(:,:,100);

z = smoothn({U_DaVis_100, V_DaVis_100});
u_smooth_100 = z{1,1}; 
v_smooth_100 = z{1,2}; 

window = 6; %for SRP images
step = 3;

n = size(data_PIV_y(:,:,100),2);%number of columns of each frame(matrix)

u = zeros(size(Matx(:,:))); %Initial u, v matrixes for each frame 
v = zeros(size(Maty(:,:)));

limit = 1;
diff_limit_lower = 0.05;
diff_limit_upper = 0.5;

for ii = 1:n
    %Define limits of w - columns of OF data
    left_s = (ii*step) - 2;
    right_s = (ii*step);
    
%     left_s = (ii-1)*step + 1;
%     right_s = (ii*step);
    
    for w = left_s:right_s
        %[u(:,w),v(:,w)] = Scaling_mean_DaVisOF(Matx(:,w),Maty(:,w), data_PIV_x(:,ii,100),data_PIV_y(:,ii,100), limit, diff_limit_lower, diff_limit_upper,window,step);
        [u(:,w),v(:,w)] = Scaling_mean_DaVisOF(Matx(:,w),Maty(:,w), u_smooth_100(:,ii),v_smooth_100(:,ii), limit, diff_limit_lower, diff_limit_upper,window,step);
    end
    
end


toc