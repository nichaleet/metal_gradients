function calibrate_flux_toHSTcat

f=load('model_colors_cleaned10052013.txt');
g=load('a383_f475wnotext.cat');

for i=1:length(f(:,1))
    for j=1:length(g(:,1))
   if abs(f(i,2)-g(j,4))<6 && abs(f(i,3)-g(j,5))<6 
       f(i,20)=g(j,6); %475 cal
       f(i,21)=g(j,74);
   end
    end
end
g=load('a383_f814wnotext.cat');
for i=1:length(f(:,1))
    for j=1:length(g(:,1))
   if abs(f(i,2)-g(j,4))<3 && abs(f(i,3)-g(j,5))<3 
       f(i,22)=g(j,6); %814 cal
       f(i,23)=g(j,74);
   end
    end
end
g=load('a383_f625wnotext.cat');
for i=1:length(f(:,1))
    for j=1:length(g(:,1))
   if abs(f(i,2)-g(j,4))<3 && abs(f(i,3)-g(j,5))<3 
       f(i,24)=g(j,6); %814 cal
       f(i,25)=g(j,74);
   end
    end
end
save 'cal_cat.txt' f -ascii