function [meanImage] = mean_image(fileNameRef,fileName, firstFrame, lastFrame,num)


% fileNameRef = 'TEST_00.tif'
% fileName = 'TEST_0'

I0 = imread(fileNameRef); % fileNameRef = 'TEST_00'
sumImage = double(I0); % Inialize to first image.
for i=firstFrame:lastFrame % Read in remaining images.
  rgbImage = imread([fileName,num2str(i),'.tif']); %fileName = 'TEST_0'
  sumImage = sumImage + double(rgbImage);
end
meanImage = sumImage /num;


