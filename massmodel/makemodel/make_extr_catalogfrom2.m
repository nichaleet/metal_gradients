function make_extr_catalogfrom2;
f=load('a383_f475wnotext.cat');
g=load('a383_f814wnotext.cat');
lines1=length(f(:,1));
for i=1:lines1
  
 new_file(i,1)=f(i,1);
 new_file(i,2)=f(i,4);
 new_file(i,3)=f(i,5);
 new_file(i,4)=f(i,6);%475
 new_file(i,5)=f(i,74);%475 10^((f(i,33)-26)/(-2.5));
  new_file(i,6)=g(i,6);%775
 new_file(i,7)=g(i,74);%775 10^((f(i,38)-26)/(-2.5));
%   new_file(i,8)=f(i,25);%775
%  new_file(i,9)=10^((f(i,25)-23)/(-2.5));
%   new_file(i,10)=f(i,30);%625
%  new_file(i,11)=10^((f(i,30)-23)/(-2.5));
 %fwhm:
 new_file(i,8)= new_file(i,4)-new_file(i,6);
new_file(i,9)=g(i,75);%-f(i,30); fwhm
new_file(i,10)=0;%f(i,8);%-f(i,30); stell
% new_file(i,11)=f(i,68);%-f(i,30); 775 mag
% new_file(i,12)=f(i,71);%-f(i,30); 775 flux
new_file(i,11)=f(i,2);%-f(i,30); ra
new_file(i,12)=f(i,3);%-f(i,30); dec
%new_file(i,10)=f(i,8);%-f(i,10); 
   end
save temp_cat_all.txt new_file -ascii;
%file_data=load('model_colors.txt');
%scatter(file_data(:,7),file_data(:,8));