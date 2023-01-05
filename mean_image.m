function [meanImage] = mean_image(filePattern, num)

%First frame is always 1 and last frame corresponds to the total number of
%frames

ImgSeq = readImgSeq(filePattern,1,num);
%num is the number of image to process
%reference image I_ref is the one with 000

I_ref = ImgSeq(:,:,1);

sumImage = double(I_ref); % Inialize to first image.
for i=2:num % Read in remaining images.
  nextImage = ImgSeq(:,:,i);
  sumImage = sumImage + double(nextImage);
end

meanImage = sumImage /num;



