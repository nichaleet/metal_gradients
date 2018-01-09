function clean_catalog_frommanual(filename);
file_data=load(filename);
g=load('GalsAmatanotextADI.reg');

%first match objects:
for i=1:length(g(:,1))
   for j=1:length(file_data(:,1))
      if sqrt((g(i,1)-file_data(j,2))^2+(g(i,2)-file_data(j,3))^2)<7
       new_file(i,:)=file_data(j,:);
       
    
      end
   end
end

% i=1;
% for obj=1:length(file_data(:,1))
%    if ((file_data(obj,8)>=1.6) && (file_data(obj,8)<2.8)) %&& (file_data(obj,6)>=15.5) && (file_data(obj,6)<=23))
%         if ((file_data(obj,6)>=19.2) && (file_data(obj,6)<23))
%         if ((file_data(obj,10)<0.9))% && (file_data(obj,8)<=0.72*file_data(obj,7)+0.1) && (file_data(obj,8)>=0.72*file_data(obj,7)-0.1))
%           if ((file_data(obj,9)>0.2)) %&&  (file_data(obj,6)>15))
%             new_file(i,:)=file_data(obj,:);
%             i=i+1;
%           end
%         end
%         end
%     end
% obj
% end

save model_colors_temp.txt new_file -ascii;