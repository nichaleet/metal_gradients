function filter_red_sequence
%enter two points and width
f=load('temp_cat_all.txt');
%scatter(f(:,6),f(:,8))
x1=12.99; y1=2.817;
x2=16.96; y2=0.8197;
widthi=1.5;

a=(y2-y1)/(x2-x1);
b=y2-a*x2;

counti=0;
for i=1:length(f(:,1))
   ytemp=f(i,6)*a+b; 
       if abs(ytemp-f(i,8))<widthi
           %cut in mag
        if ((f(i,6)>=10) && (f(i,6)<18)) % magnitude limit
            %stellar obj?
         if ((f(i,9)>0.8)) %&&  (file_data(obj,6)>15))
           % if  f(i,15)==-99 %&& (abs(f(i,15)-0.308)<0.06) || (file_data(obj,6)>15))
           
             counti=counti+1;
           new_file(counti,:)=f(i,:);
           indi(counti)=i;
            end
          %end
        end
       end
end
scatter(f(:,6),f(:,8)); hold on;
xvec=[x1:0.01:x2];
yvec=a.*xvec+b; 
plot(xvec,yvec,'black');
plot(xvec,yvec+widthi,'black');
plot(xvec,yvec-widthi,'black');
display('click on image to finish and see chosen members')
waitforbuttonpress
scatter(f(indi,6),f(indi,8),'red');
hold off
save model_colors_templim24.txt new_file -ascii