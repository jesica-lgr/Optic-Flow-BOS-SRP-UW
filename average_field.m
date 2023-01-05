function [av_field_x,av_field_y] = average_field(data_x, data_y, center,num_frames)
%Obtains an averaged field considering the two subsequent and prior frames.
%num_frames must be an odd number

[m,n,z] = size(data_x);

av_field_x = zeros(m,n);
av_field_y = zeros(m,n);

for i = 1:m
    for j = 1:n
        num= floor(num_frames/2);
        sum_data_x = 0;
        sum_data_y = 0;
        for k =(center-num):(center+num)
            sum_data_x = sum_data_x + data_x(i,j,k);
            sum_data_y = sum_data_y + data_y(i,j,k);
        end
        
        av_field_x(i,j) = sum_data_x/num_frames;
        av_field_y(i,j) = sum_data_y/num_frames;
%         av_field_x(i,j) = [data_x(i,j,98)+ data_x(i,j,99)+data_x(i,j,100) + data_x(i,j,101)+data_x(i,j,102)]/5;
%         av_field_y(i,j) = [data_y(i,j,98)+ data_y(i,j,99)+data_y(i,j,100) + data_y(i,j,101)+data_y(i,j,102)]/5;
    end
end

