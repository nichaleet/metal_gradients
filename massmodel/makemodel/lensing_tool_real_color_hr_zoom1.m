function lensing_tool_real_color_hr_zoom1(k_dlds,xc,yc,sizec);
% this program will show the image, let you mark the image you are
% interested in and then will de-lens and re-lens it. "k" is proportional to Dls/Ds. 
%Image should be constructed from the original data images.
%!!!! MAGNIFICATION MUST BE CALCULATED WITH THE RIGHT FACTORS BEFOER!!!!!!
%!!!you should run "transform image" before on the image. Note that this changes with different images.
display('loading...')
load transformed10000_0647CovTRY2_2

clear ZIX ZIY alpha_x alpha_y alpha_x_DM alpha_y_DM file_data magnification magnification_DM magnification_im verification verification_DM
k=k_dlds;
display('trimming image to right size of deflection field.')
display('soon you will see the image, mark the region with mouse. Double click to finish')
%if you don't know the right transform:
x_start=i_x;
y_start=i_y;
xc=xc-x_start;
yc=yc-y_start;
try_im=new_im(xc-round(sizec/2):xc+round(sizec/2),yc-round(sizec/2):yc+round(sizec/2),:);


BW=roipolyold(try_im); 
display('OK')

BWnew(1:3000,1:3000,:)=0;
BWnew(xc-round(sizec/2):xc+round(sizec/2),yc-round(sizec/2):yc+round(sizec/2),:)=BW(:,:,:);

s_check_realm(1:3000,1:3000,1:3)=0;
source_im=new_im;
toshow=new_im;
source_im(:,:,:)=0;
toshow(:,:,:)=0;


%now activating image to source
display('de-lensing...')
for m=1:x_size-i_x-5
 for n=1:y_size-5-i_y
   if (BWnew(m,n)>0)
      
    im_x=round((m-alpha_x_ALL(m,n)*k));
    im_y=round((n-alpha_y_ALL(m,n)*k));  
    %if (new_im(m,n,:)>55)
    s_check_realm(im_x,im_y,:)=new_im(m,n,:); 
    source_im(im_x,im_y,:)=new_im(m,n,:); 
    %marking the source position
    new_im(round(im_x),round(im_y),:)=255;
%    if (new_im(m,n,:)>50)
%     round((im_x+4*750)/4)
%    round((im_y+4*750)/4)
    %end
  end
 end
end


%re-lens
median(median(new_im(:,:,1)))
median(median(new_im(:,:,2)))
median(median(new_im(:,:,3)))

display('relensing...')
for m=1:x_size-i_x-5
   for n=1:y_size-5-i_y
    im_x=round((m-alpha_x_ALL(m,n)*k));
    im_y=round((n-alpha_y_ALL(m,n)*k));  
    %new_im(m,n,:)=s_check_realm(im_x,im_y,:);
   if (im_x<6000 && im_y< 6000 && im_x>0 && im_y>0)
    %if (s_check_realm(im_x,im_y,:)>55)
    %new_im(m,n,:)= s_check_realm(im_x,im_y,:);
%     m
%     n
    new_im(m,n,1)=s_check_realm(im_x,im_y,1);
    new_im(m,n,2)=s_check_realm(im_x,im_y,2);
    new_im(m,n,3)=s_check_realm(im_x,im_y,3);
    toshow(m,n,1)=s_check_realm(im_x,im_y,1);
    toshow(m,n,2)=s_check_realm(im_x,im_y,2);
    toshow(m,n,3)=s_check_realm(im_x,im_y,3);
    %end
    end
    
    end
end
%save tempo;
% size(new_im)
% size(s_check_realm)

imwrite(source_im,'source.tif');
imwrite(new_im,'to_art.tif');
imwrite(toshow,'toshow.tif');
imshow(toshow); 
set(gca,'XDir','normal');
