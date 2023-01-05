function [data_x, data_y] = Scaling_mean(data_Horn_x,data_Horn_y, data_PIV_x,data_PIV_y, limit, diff_limit_lower, diff_limit_upper, window,step)

%limit refers to the PIV value from which we assumme that  scaling is needed in principle
%diff_limit refres to the max difference between Horn's and PIV data 
%This function results focuses on the begining and end of the sets and
%appliesa tree point mean to smooth the data
%Al data sets are assumed to be of size mx1


data_interval = 0.1;
width = 8;
[m,n] = size(data_PIV_x);


data_x = zeros(size(data_Horn_x,1)-(step),1);
data_y = zeros(size(data_Horn_y,1)-(step),1);


%Get general scale
scale_x = Myscale1(data_Horn_x, data_PIV_x,data_interval);
scale_y = Myscale1(data_Horn_y, data_PIV_y,data_interval);


%k corresponds to the value in data_PIV
for k = 1:m
    
    %Define range of Horn's set 
    %left_val = (k*8)-7;
    %right_val = k*8;
    clear left_val_1;
    clear left_val_2;
    clear right_val_1;
    clear right_val_1;
    
    %left_val_1 = (k*8)-7;%+4
    %right_val_2 = k*8;
    
    %right_val_1 = (k*8)/2;
    %left_val_2 = ((k*8)/2)+1;
    %right_val_1 = left_val_1+3;
    %left_val_2 = right_val_2-3;
    
%     left_val_1 = (k*8)-3;
%     right_val_2 = k*8+4;
%     right_val_1 = left_val_1+3;
%     left_val_2 = right_val_2-3;
    
    left_val_1 = (k*step)-((step/2)-1);
    right_val_2 = k*step+(step/2);
    %right_val_2 = k*step;
    right_val_1 = left_val_1+((step/2)-1);
    left_val_2 = right_val_2-((step/2)-1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%X VALUES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if isnan(data_PIV_x(k,1))% Llenar ese espacio con el promedio de los dos vecinos 
        if k == 1 
            if isnan(data_PIV_x(k+1,1))
                data_PIV_x(k,1) = data_PIV_x(k+2,1);
            else
                data_PIV_x(k,1) = data_PIV_x(k+1,1);
            end
        elseif k == m
            data_PIV_x(k,1) = data_PIV_x(k-1,1); 
        elseif k== m-1
            if isnan(data_PIV_x(k+1,1)) && ~isnan(data_PIV_x(k-1,1))
                data_PIV_x(k,1) = data_PIV_x(k-1,1);
            elseif isnan(data_PIV_x(k-1,1)) && ~isnan(data_PIV_x(k+1,1))
                data_PIV_x(k,1) = (data_PIV_x(k-2,1)+ data_PIV_x(k+1,1))/2;
            else
                data_PIV_x(k,1) = (data_PIV_x(k-1,1)+ data_PIV_x(k+1,1))/2;
            end
        elseif k== m-2
            if isnan(data_PIV_x(k+1,1)) && ~isnan(data_PIV_x(k-1,1))
                data_PIV_x(k,1) = data_PIV_x(k-1,1);
            elseif isnan(data_PIV_x(k-1,1)) && ~isnan(data_PIV_x(k+1,1))
                data_PIV_x(k,1) = (data_PIV_x(k-2,1)+ data_PIV_x(k+1,1))/2;
            else
                data_PIV_x(k,1) = (data_PIV_x(k-1,1)+ data_PIV_x(k+1,1))/2;
            end
        else
            if isnan(data_PIV_x(k+1,1)) && ~isnan(data_PIV_x(k+2,1))
                data_PIV_x(k,1) = (data_PIV_x(k-1,1)+ data_PIV_x(k+2,1))/2;
                %data_PIV_x(k,1) = data_PIV_x(k-1,1);
            elseif isnan(data_PIV_x(k-1,1)) && ~isnan(data_PIV_x(k+1,1))
                data_PIV_x(k,1) = (data_PIV_x(k-2,1)+ data_PIV_x(k+1,1))/2;
            else
                data_PIV_x(k,1) = (data_PIV_x(k-1,1)+ data_PIV_x(k+1,1))/2;
            end
        end
    end
        
    
    if abs(data_PIV_x(k,1)) <= limit %In principle, no scaling is needed
        
        
        for i = left_val_1:right_val_1 %Extra condition considering Horn's set of 8 values
            diff_Horn_PIV_x = abs(data_PIV_x(k,1) - data_Horn_x(i,1)); %Difference between Horn and PIV
            
            if diff_Horn_PIV_x > diff_limit_lower %Scaling is needed
                n_interval_x = data_Horn_x(i,1)/data_interval;
                data_x(i,1) = data_Horn_x(i,1) + n_interval_x*scale_x*0.4;%Magnification
                
                %Sum or sustract the diference
                if data_PIV_x(k,1) < 0
                    data_x(i,1) = data_x(i,1) - diff_Horn_PIV_x;
                else
                    data_x(i,1) = data_x(i,1) + diff_Horn_PIV_x;
                end     
                
            else %No scaling is needed
                data_x(i,1) = data_Horn_x(i,1);
            end 
        end
        
        for i = left_val_2:right_val_2 %Extra condition considering Horn's set of 8 values
            diff_Horn_PIV_x = abs(data_PIV_x(k,1) - data_Horn_x(i,1)); %Difference between Horn and PIV
            
            if diff_Horn_PIV_x > diff_limit_lower %Scaling is needed
                n_interval_x = data_Horn_x(i,1)/data_interval;
                data_x(i,1) = data_Horn_x(i,1) + n_interval_x*scale_x*0.4;%Magnification
                
                %Sum or sustract the diference
                if data_PIV_x(k,1) < 0
                    data_x(i,1) = data_x(i,1) - diff_Horn_PIV_x;
                else
                    data_x(i,1) = data_x(i,1) + diff_Horn_PIV_x;
                end     
                
            else %No scaling is needed
                data_x(i,1) = data_Horn_x(i,1);
            end 
          
        end
         
            
 %I am assuming that if abs(data_PIV_x(k,1)) > limit, automatically scaling
 %is required ???
        
    elseif abs(data_PIV_x(k,1)) > limit %Scaling is needed  
        
        
        for i = left_val_1:right_val_1 %Extra condition considering Horn's set of 8 values
            diff_Horn_PIV_x = abs(data_PIV_x(k,1) - data_Horn_x(i,1)); %Difference between Horn and PIV
            
            if diff_Horn_PIV_x > diff_limit_upper %Scaling is needed
                n_interval_x = data_Horn_x(i,1)/data_interval;
                data_x(i,1) = data_Horn_x(i,1) + n_interval_x*scale_x*0.4;%Magnification
                
                %Sum or sustract the diference
                if data_PIV_x(k,1) < 0
                    data_x(i,1) = data_x(i,1) - diff_Horn_PIV_x;
                else
                    data_x(i,1) = data_x(i,1) + diff_Horn_PIV_x;
                end     
                
            else %No scaling is needed
                data_x(i,1) = data_Horn_x(i,1);
            end 
        end  
        
        for i = left_val_2:right_val_2 %Extra condition considering Horn's set of 8 values
            diff_Horn_PIV_x = abs(data_PIV_x(k,1) - data_Horn_x(i,1)); %Difference between Horn and PIV
            
            if diff_Horn_PIV_x > diff_limit_lower %Scaling is needed
                n_interval_x = data_Horn_x(i,1)/data_interval;
                data_x(i,1) = data_Horn_x(i,1) + n_interval_x*scale_x*0.4;%Magnification
                
                %Sum or sustract the diference
                if data_PIV_x(k,1) < 0
                    data_x(i,1) = data_x(i,1) - diff_Horn_PIV_x;
                else
                    data_x(i,1) = data_x(i,1) + diff_Horn_PIV_x;
                end     
                
            else %No scaling is needed
                data_x(i,1) = data_Horn_x(i,1);
            end 
        end
           
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%Y VALUES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isnan(data_PIV_y(k,1))% Llenar ese espacio con el promedio de los dos vecinos 
        if k == 1 
            if isnan(data_PIV_y(k+1,1))
                data_PIV_y(k,1) = data_PIV_y(k+2,1);
            else
                data_PIV_y(k,1) = data_PIV_y(k+1,1);
            end
        elseif k == m
            data_PIV_y(k,1) = data_PIV_y(k-1,1);
           
        elseif k== m-1
            if isnan(data_PIV_y(k+1,1)) && ~isnan(data_PIV_y(k-1,1))
                data_PIV_y(k,1) = data_PIV_y(k-1,1);
            elseif isnan(data_PIV_y(k-1,1)) && ~isnan(data_PIV_y(k+1,1))
                data_PIV_y(k,1) = (data_PIV_y(k-2,1)+ data_PIV_y(k+1,1))/2;
            else
                data_PIV_y(k,1) = (data_PIV_y(k-1,1)+ data_PIV_y(k+1,1))/2;
            end
        elseif k== m-2
            if isnan(data_PIV_y(k+1,1)) && ~isnan(data_PIV_y(k-1,1))
                data_PIV_y(k,1) = data_PIV_y(k-1,1);
            elseif isnan(data_PIV_y(k-1,1)) && ~isnan(data_PIV_y(k+1,1))
                data_PIV_y(k,1) = (data_PIV_y(k-2,1)+ data_PIV_y(k+1,1))/2;
            else
                data_PIV_y(k,1) = (data_PIV_y(k-1,1)+ data_PIV_y(k+1,1))/2;
            end
        else
            if isnan(data_PIV_y(k+1,1)) && ~isnan(data_PIV_y(k+2,1))
                data_PIV_y(k,1) = (data_PIV_y(k-1,1)+ data_PIV_y(k+2,1))/2;
                %data_PIV_x(k,1) = data_PIV_x(k-1,1);
            elseif isnan(data_PIV_y(k-1,1)) && ~isnan(data_PIV_y(k+1,1))
                data_PIV_y(k,1) = (data_PIV_y(k-2,1)+ data_PIV_y(k+1,1))/2;
            else
                data_PIV_y(k,1) = (data_PIV_y(k-1,1)+ data_PIV_y(k+1,1))/2;
            end
        end
    end
      
    
    if abs(data_PIV_y(k,1)) <= limit %In principle, no scaling is needed 
        
        for i = left_val_1:right_val_1 %Extra condition considering Horn's set of 8 values
            diff_Horn_PIV_y = abs(data_PIV_y(k,1) - data_Horn_y(i,1)); %Difference between Horn and PIV
            
            if diff_Horn_PIV_y > diff_limit_lower %Scaling is needed
                n_interval_y = data_Horn_y(i,1)/data_interval;
                data_y(i,1) = data_Horn_y(i,1) + n_interval_y*scale_y*0.4;%Magnification
                
                %Sum or sustract the diference
                if data_PIV_y(k,1) < 0
                    data_y(i,1) = data_y(i,1) - diff_Horn_PIV_y;
                else
                    data_y(i,1) = data_y(i,1) + diff_Horn_PIV_y;
                end     
                
            else %No scaling is needed
                data_y(i,1) = data_Horn_y(i,1);
            end   
        end
        
        for i = left_val_2:right_val_2 %Extra condition considering Horn's set of 8 values
            diff_Horn_PIV_y = abs(data_PIV_y(k,1) - data_Horn_y(i,1)); %Difference between Horn and PIV
            
            if diff_Horn_PIV_y > diff_limit_lower %Scaling is needed
                n_interval_y = data_Horn_y(i,1)/data_interval;
                data_y(i,1) = data_Horn_y(i,1) + n_interval_y*scale_y*0.4;%Magnification
                
                %Sum or sustract the diference
                if data_PIV_y(k,1) < 0
                    data_y(i,1) = data_y(i,1) - diff_Horn_PIV_y;
                else
                    data_y(i,1) = data_y(i,1) + diff_Horn_PIV_y;
                end     
                
            else %No scaling is needed
                data_y(i,1) = data_Horn_y(i,1);
            end 
        end
  

 %I am assuming that if abs(data_PIV_y(k,1)) > limit, automatically scaling
 %is required ???
        
    elseif abs(data_PIV_y(k,1)) > limit %Scaling is needed  
        
        for i = left_val_1:right_val_1 %Extra condition considering Horn's set of 8 values
            diff_Horn_PIV_y = abs(data_PIV_y(k,1) - data_Horn_y(i,1)); %Difference between Horn and PIV
            
            if diff_Horn_PIV_y > diff_limit_upper %Scaling is needed
                n_interval_y = data_Horn_y(i,1)/data_interval;
                data_y(i,1) = data_Horn_y(i,1) + n_interval_y*scale_y*0.4;%Magnification
                
                %Sum or sustract the diference
                if data_PIV_y(k,1) < 0
                    data_y(i,1) = data_y(i,1) - diff_Horn_PIV_y;
                else
                    data_y(i,1) = data_y(i,1) + diff_Horn_PIV_y;
                end     
                
            else %No scaling is needed
                data_y(i,1) = data_Horn_y(i,1);
            end 
        end
        
        for i = left_val_2:right_val_2 %Extra condition considering Horn's set of 8 values
            diff_Horn_PIV_y = abs(data_PIV_y(k,1) - data_Horn_y(i,1)); %Difference between Horn and PIV
            
            if diff_Horn_PIV_y > diff_limit_lower %Scaling is needed
                n_interval_y = data_Horn_y(i,1)/data_interval;
                data_y(i,1) = data_Horn_y(i,1) + n_interval_y*scale_y*0.4;%Magnification
                
                %Sum or sustract the diference
                if data_PIV_y(k,1) < 0
                    data_y(i,1) = data_y(i,1) - diff_Horn_PIV_y;
                else
                    data_y(i,1) = data_y(i,1) + diff_Horn_PIV_y;
                end     
                
            else %No scaling is needed
                data_y(i,1) = data_Horn_y(i,1);
            end 
        end 
        
    end  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Work with the next set of 8 values
    k = k+1;   
    
end


%fix begining and end
        
datax_mean = data_x;
datay_mean = data_y;

% for k_mean = 1:m-1
for k_mean = 1:m-1
    
%     left_val = (k_mean*8)-3;
%     right_val = k_mean*8+4;

      left_val = (k_mean*step)-((step/2)-1);
      right_val = k_mean*step+(step/2);
% 
%     left_val = (k_mean*step)-(step-1);
%     right_val = k_mean*step;
    
    datax_mean(left_val,1) = (data_x(left_val-1,1) + data_x(left_val,1) + data_x(left_val+1,1))/3;
    datax_mean(right_val,1) = (data_x(right_val-1,1) + data_x(right_val,1) + data_x(right_val+1,1))/3;
    
    datay_mean(left_val,1) = (data_y(left_val-1,1) + data_y(left_val,1) + data_y(left_val+1,1))/3;
    datay_mean(right_val,1) = (data_y(right_val-1,1) + data_y(right_val,1) + data_y(right_val+1,1))/3;
    
    %datax_mean(left_val,1) = (datax_mean(left_val-1,1) + datax_mean(left_val,1) + datax_mean(left_val+1,1))/3;
    %datax_mean(right_val,1) = (datax_mean(right_val-1,1) + datax_mean(right_val,1) + datax_mean(right_val+1,1))/3;
    
    %datay_mean(left_val,1) = (datay_mean(left_val-1,1) + datay_mean(left_val,1) + datay_mean(left_val+1,1))/3;
    %datay_mean(right_val,1) = (datay_mean(right_val-1,1) + datay_mean(right_val,1) + datay_mean(right_val+1,1))/3;
    
end










%Normalize data 

%data_x_norm = scale_adjust_minmax(data_x);
%data_y_norm = scale_adjust_minmax(data_y);
%data_PIV_xnorm = scale_adjust_minmax(data_PIV_x);
%data_PIV_ynorm = scale_adjust_minmax(data_PIV_y);


%Convert to meters and plot data_Horn_scaled and data_PIV

% figure(1)
% plot_meters_scatter(data_x,1)
% %plot_meters(data_x,1)
% hold on; 
% plot_meters_scatter(data_PIV_x,8); 
% hold off
% 
% figure(2)
% plot_meters_scatter(data_y,1)
% %plot_meters(data_y,1)
% hold on; 
% plot_meters_scatter(data_PIV_y,8); 
% hold off
% 
% figure(3)
% plot_meters_scatter(datax_mean,1)
% %plot_meters(datax_mean,1)
% hold on; 
% plot_meters_scatter(data_PIV_x,8); 
% hold off
% 
% figure(4)
% plot_meters_scatter(datay_mean,1)
% %plot_meters(datay_mean,1)
% hold on; 
% plot_meters_scatter(data_PIV_y,8); 
% hold off








%figure(3)

%plot_meters_scatter(data_x_norm,1)
%plot_meters(data_x_norm,1)
%hold on; 
%plot_meters_scatter(data_PIV_xnorm,8); 
%hold off


%figure(4)
%plot_meters_scatter(data_y_norm,1)
%plot_meters(data_y_norm,1)
%hold on; 
%plot_meters_scatter(data_PIV_ynorm,8); 
%hold off
end