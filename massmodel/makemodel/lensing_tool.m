function lensing_tool(k_dlds);
% this program will show the image, let you mark the image you are
% interested in and then will de-lens and re-lens it. "k" is proportional to Dls/Ds. 
%Image should be constructed from the original data images.
%!!!! MAGNIFICATION MUST BE CALCULATED WITH THE RIGHT FACTORS BEFOER!!!!!!
%!!!you should run "transform image" before on the image. Note that this changes with different images.
display('loading...')
%load transformed10000_NEWFILE;
%load transformed2800011495FINAL
%load transformed5200_11495photz1
%load transformed35000_11495_10extras
%load transformedignore
%load transformedinitial
load transformedSecond
%load transformedhighergal
%load transformedTry1BCG1RedApr20Gaussiancont3
%load transformedspeedup4free20000TRY4
%load transformedspeedup3free30000fixed1
%load transformed10000_4SS1good
%load transformedspeedup3free30000fixed1
clear ZIX ZIY alpha_x alpha_y alpha_x_DM alpha_y_DM file_data magnification magnification_DM magnification_im verification verification_DM
k=k_dlds;
display('trimming image to right size of deflection field.')
display('soon you will see the image, mark the region with mouse. Double click to finish')
%if you don't know the right transform:

% load magnif_ALL
% img=imread(image_name);
% 
% temp_im=img(4488-3779:4488-1278,4368-3416:4368-915,:);
% 
% for m=1:x_size-i_x
%  for n=1:y_size-i_y
%    new_im(m,n,:)=temp_im(2502-n,m,:);
%      
%  end
% end


BW=roipolyold(new_im); 
display('OK')
%mask=poly2mask(BW);

%save lensed;
imgo(1:leng-5,1:leng-5)=0;
s_check_real(1:leng-5,1:leng-5)=0;

display('copying source')

%now activating image to source
display('de-lensing...')
for m=1:leng-5
 for n=1:leng-5
   if (BW(m,n)>0)
    m
    n
    im_x=m-round(alpha_x_ALL(m,n)*k);
    im_y=n-round(alpha_y_ALL(m,n)*k);  
    s_check_real(im_x,im_y)=new_im(m,n); 
    %marking the source position
    new_im(im_x,im_y,:)=255;
   
  end
 end
end


%re-lens
display('relensing...')
for m=1:leng-5
   for n=1:leng-5
    im_x=m-round(alpha_x_ALL(m,n)*k);
    im_y=n-round(alpha_y_ALL(m,n)*k);
    if (im_x >0 && im_x<leng-10 && im_y>0 && im_y<leng-10)
    imgo(m,n)=s_check_real(im_x,im_y);
    end
    
    if (imgo(m,n)>0)
    new_im(m,n,1)=255;
    new_im(m,n,2)=0;
    new_im(m,n,3)=0;
   

   end
    
    end
end
%save tempo;

imshow(new_im); 
set(gca,'XDir','normal');
