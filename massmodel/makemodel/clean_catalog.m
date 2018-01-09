function clean_catalog(filename);
file_data=load(filename);
i=1;
for obj=1:length(file_data(:,1))
   if ((file_data(obj,8)>=1) && (file_data(obj,8)<3)) %&& (file_data(obj,6)>=15.5) && (file_data(obj,6)<=23))
        if ((file_data(obj,6)>=18) && (file_data(obj,6)<25))
        %if ((file_data(obj,12)>4))% && (file_data(obj,8)<=0.72*file_data(obj,7)+0.1) && (file_data(obj,8)>=0.72*file_data(obj,7)-0.1))
          if ((file_data(obj,10)<0.2)) %&&  (file_data(obj,6)>15))
            new_file(i,:)=file_data(obj,:);
            i=i+1;
         % end
        end
        end
    end
obj
end

save model_colors_temp.txt new_file -ascii;
%file_data=load('model_colors.txt');
%scatter(file_data(:,7),file_data(:,8));