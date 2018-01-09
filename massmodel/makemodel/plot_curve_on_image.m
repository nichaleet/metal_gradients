function plot_curve_on_image(im);
%load bestmodel_aCovtry
%load transformed10000_NEWFILE;
%load transformed2800011495FINAL
%load transformed5200_11495photz1
%load transformed35000_11495_10extras
%load transformedignore
%load transformed16052012prelim
load transformedLTMv18102013
%load transformedhighergal
%load transformedTry1BCG1RedApr20Gaussiancont3
%load transformedspeedup4free20000TRY4
 %load model_example_gaussion10000;
 %load BestModel7Free4RerdsFree20000TRY3
 %load transformedspeedup&free10000TRY3
 %load transformedspeedup3free20000fixed1
 %load transformed10000_4SS1good
 %load transformedspeedup3free30000fixed1
 %load transformedRegProc10000_3
 %load magnif_ALL_z3
 %load magnif_ALL_z2p474free20000TRY3
% 
% x_size
% i_x
% %create the curve only if > from 4 for example
% to_ring=log((magnification_im_ALL+1));
%load magnif_ALL_z1p89_v3FINAL
to_ring=log((magnification_ALL'+1));
new_im=imread(im);
for m=1:(x_size-i_x-5)
   % m
 for  n=1:(y_size-i_y-5)
     
     if (to_ring(m,n)>=6)
         
         new_im(leng-m,n,:)=255;
     end
%       if (to_ring2(m,n)>=8)
%          
%          new_im(2501-m,n,1)=0;
%          new_im(2501-m,n,2)=100;
%          new_im(2501-m,n,3)=100;
%      end    
%      
 end
end

%imwrite(new_im,'magOnIm.tif');
imshow(new_im(1:leng,1:leng,:)); 
set(gca,'XDir','normal');
%f=imread('z2adi.tif');
%size(f)
%set(gca,'YDir','normal');
%  imagesc(new_im);
% %colormap(gray)
% %colormap(gray);
% set(gca,'YDir','normal'); 