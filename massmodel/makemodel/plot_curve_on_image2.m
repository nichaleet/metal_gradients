function plot_curve_on_image2(im);
xoff=155
yoff=215
 load magnif_ALL_z3p5;
% to_ring=log((magnification_im_ALL+1));
 new_im00=imread(im);
x_size=3000
i_x=2
%create the curve only if > from 4 for example
to_ring=log((magnification_ALL'+1));
for m=1+i_x:(x_size-i_x-5)
   % m
 for  n=1+i_x:(x_size-i_x-5)
     
     if (to_ring(3000-m,n)>=6.5 && m-xoff<3000 && m-xoff>0 && n-yoff<3000 && n-yoff>0)
         new_im00(m-xoff,n-yoff,1)=200;
         new_im00(m-xoff,n-yoff,2)=50;
         new_im00(m-xoff,n-yoff,3)=50;
     end
 end
end
load transformed10000_0647CovTRY3_1
to_ring=log((magnification_ALL'+1));
x_size=3000
i_x=2
for m=1+i_x:(x_size-i_x-5)
   % m
 for  n=1+i_x:(x_size-i_x-5)
     
     if (to_ring(3000-m,n)>=6.5 && m-xoff<3000 && m-xoff>0 && n-yoff<3000 && n-yoff>0)
         
         new_im00(m-xoff,n-yoff,:)=255;
     end
 end
end
load magnif_ALL_z5
to_ring=log((magnification_ALL'+1));
x_size=3000
i_x=2
for m=1+i_x:(x_size-i_x-5)
   % m
 for  n=1+i_x:(x_size-i_x-5)
     
     if (to_ring(3000-m,n)>=6.5 && m-xoff<3000 && m-xoff>0 && n-yoff<3000 && n-yoff>0)
         
         new_im00(m-xoff,n-yoff,1)=50;
         new_im00(m-xoff,n-yoff,2)=50;
         new_im00(m-xoff,n-yoff,3)=200;
         
     end
 end
end
load magnif_ALL_z11
to_ring=log((magnification_ALL'+1));
x_size=3000
i_x=2
for m=1+i_x:(x_size-i_x-5)
   % m
 for  n=1+i_x:(x_size-i_x-5)
     
     if (to_ring(3000-m,n)>=6.5 && m-xoff<3000 && m-xoff>0 && n-yoff<3000 && n-yoff>0)
         
         new_im00(m-xoff,n-yoff,1)=255;
         new_im00(m-xoff,n-yoff,2)=0;
         new_im00(m-xoff,n-yoff,3)=0;
         
     end
 end
end
% % load magnif_ALL_z8
% % to_ring=log((magnification_im_ALL+1));
% % x_size=3000
% % i_x=200
% % for m=1+i_x:(x_size-i_x-5)
% %    % m
% %  for  n=1+i_x:(x_size-i_x-5)
% %      
% %      if (to_ring(3000-m,n)>=6.7)
% %          
% %          new_im00(m,n,1)=255;
% %          new_im00(m,n,2)=0;
% %          new_im00(m,n,3)=0;
% %      end
% %  end
% % end

imwrite(new_im00,'4curvesOnIm0647vMC.tif');
imshow(new_im00(:,:,:)); 
set(gca,'XDir','normal');
%f=imread('z2adi.tif');
%size(f)
%set(gca,'YDir','normal');
%  imagesc(new_im);
% %colormap(gray)
% %colormap(gray);
% set(gca,'YDir','normal'); 