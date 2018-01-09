function plot_circles_like_regions(im,model);

i_x=1250;
i_y=1250;
img=imread(im);
b=size(img);
f=load(model);
b
%img2=img;
%img2(1:b(1),1:b(2))=0;
l=length(f(:,1))
 for obj=1:l
     if ((round(f(obj,2)-i_x) >10 ) &&(round(f(obj,2)-i_x) <=2500-5 ) &&  (round(2501-f(obj,3)+i_y)>10 )&&(round(2501-f(obj,3)+i_y) <= 2500-5)) 
        %f(obj,3) 
       %img(round(f(obj,2)-i_x)-5:round(f(obj,2)-i_x)+5,round(2501-f(obj,3)+i_y)-5:round(2501-f(obj,3)+i_y)+5,1)=255;
       img(round(2501-f(obj,3)+i_y)-5:round(2501-f(obj,3)+i_y)+5, round(f(obj,2)-i_x)-5:round(f(obj,2)-i_x)+5,1)=255;
       img(round(2501-f(obj,3)+i_y)-5:round(2501-f(obj,3)+i_y)+5, round(f(obj,2)-i_x)-5:round(f(obj,2)-i_x)+5,2)=0;
       img(round(2501-f(obj,3)+i_y)-5:round(2501-f(obj,3)+i_y)+5, round(f(obj,2)-i_x)-5:round(f(obj,2)-i_x)+5,3)=0;

     end        
 end
% mask=img2>0;
  %size(img2)
%out1 = imoverlay(img, mask, [1 0 0]);
imagesc(img);
