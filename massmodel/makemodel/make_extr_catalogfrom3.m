function make_extr_catalogfrom3;
f=load('475byAdinotext.cat');
g=load('814byAdinotext.cat');
h=load('625byAdinotext.cat');
lines1=length(f(:,1));
for i=1:1101:315:1416
  
 new_file(i,1)=f(i,1);
 new_file(i,2)=f(i,2);
 new_file(i,3)=f(i,3);
 new_file(i,4)=f(i,6);%475
 new_file(i,5)=f(i,12);%475 10^((f(i,33)-26)/(-2.5));
  new_file(i,6)=g(i,6);%814
 new_file(i,7)=g(i,12);%814 10^((f(i,38)-26)/(-2.5));
  new_file(i,8)=h(i,6);%625
  new_file(i,9)=h(i,12);%625
%   new_file(i,10)=f(i,30);%625
%  new_file(i,11)=10^((f(i,30)-23)/(-2.5));
 %fwhm:
 new_file(i,10)= new_file(i,4)-new_file(i,6);
  new_file(i,11)= new_file(i,8)-new_file(i,6);
new_file(i,12)=g(i,18);%-f(i,30); fwhm
new_file(i,13)=g(i,25);%-f(i,30); star/galaxy flag
% new_file(i,11)=f(i,68);%-f(i,30); 775 mag
% new_file(i,12)=f(i,71);%-f(i,30); 775 flux
new_file(i,14)=f(i,4);%-f(i,30); ra
new_file(i,15)=f(i,5);%-f(i,30); dec
new_file(i,16:19)=f(i,21:24);%-f(i,30); ellipticity params: a,b,theta,elli

%new_file(i,10)=f(i,8);%-f(i,10); 
end
   %new_file
%save temp_cat_all.txt new_file -ascii;
save ignore.txt new_file -ascii;
%file_data=load('model_colors.txt');
%scatter(file_data(:,7),file_data(:,8));