function scale = Myscale1(data_Horn, data_PIV,data_interval)

%data_Horn and data_PIV must be vectors of mx1


[m,n]=size(data_Horn);

interval_Horn = max(data_Horn)-min(data_Horn);
interval_PIV = max(data_PIV)-min(data_PIV);



scale = data_interval*(interval_PIV/interval_Horn); %for each data_interval in Horn, there is an increment
%of x in the PIVdata