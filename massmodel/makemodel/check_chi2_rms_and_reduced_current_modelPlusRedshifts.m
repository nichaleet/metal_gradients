function check_chi2_rms_and_reduced_current_modelPlusRedshifts(pix_scale,images_file2,z_l,z_s_main)

% here need to enter only one system file as input!! no more!
%load bestmodel_aCovtry
load transformed0201
rounds=10;
minus=0.08;
final_chi2_pixels=0;
rms=0;
dg=0.01;%minus*2/rounds;
for ing=1:rounds
images_file=images_file2;
display('delensing-relensing and calculating chi^2')
%now calculatechi2
%[final_chi2_pixels]=lensing_toolBestSim3chirxj1347('allimages1206nocandidates.txt');
%[final_chi2_pixels]=lensing_toolBestSim3chirxj1347(images_file);
f=load(images_file);
sigma_location=0.5/pix_scale;
min_disto(1:length(f(:,1)))=10000;
count=0;
flag_lensingok=1;
count=count+1;
no_systems=f(1,5);
line_no=1;
%sum_disto(1:no_systems)=0;
for i=1:no_systems

% now working on first system
no_images_in_sys=f(line_no,4);
s_check(1:x_size-i_x,1:y_size-i_y)=0;
countr=0;
for no_images=1:no_images_in_sys
  if flag_lensingok==1  
    countr=countr+1;
clear imo 
imo(1:x_size+1-i_x,1:y_size+1-i_y)=0;
a=f(line_no,1)-x_start+1;
b=f(line_no,2)-y_start+1;

    imo(a-1:a+1,b-1:b+1)=100;
 
for m=a-3:a+3
 for n=b-3:b+3
     im_x(countr)=m-round(alpha_x_ALL(m,n)*(f(line_no,3)-minus+dg*(ing-1)));
     im_y(countr)=n-round(alpha_y_ALL(m,n)*(f(line_no,3)-minus+dg*(ing-1))) ;
     if im_x(countr)>0 && im_x(countr)<(x_size-i_x-2) && im_y(countr)>0 && im_y(countr)<(y_size-i_y-2)
     s_check(im_x,im_y)=imo(m,n); 
    s_check(im_x(countr),im_y(countr))=imo(m,n); 
         else
         display('could not lens source - out of FOV; with:')
         q,s,k_new,k_gal,gamma,phi
         flag_lensingok=0;
     end
 end
end
 end
    line_no=line_no+1;
  
end

%%until constructed source image for system now average
line_no=line_no-1;
if flag_lensingok==1
%line_no
imxavg=mean(im_x);
imyavg=mean(im_y);
s_check(1:x_size-i_x,1:y_size-i_y)=0;
s_check(round(imxavg)-1:round(imxavg)+1,round(imyavg)-1:round(imyavg)+1)=100;
clear im_x im_y
%here should measure the distance from real image
%for ja=1:no_images_in_sys
im_x=1;
im_y=1;
for m=1:leng-5
   for n=1:leng-5
%        if m==1 && n==1
%            display('started imageplane')
%        end
   
    im_x=m-round(alpha_x_ALL(m,n)*(f(line_no,3)-minus+dg*(ing-1)));
    im_y=n-round(alpha_y_ALL(m,n)*(f(line_no,3)-minus+dg*(ing-1)));
    if (im_x >0 && im_x<leng-10 && im_y>0 && im_y<leng-10)
        if s_check(im_x,im_y)>0
         
         imo(m,n)=s_check(im_x,im_y);
     
         for ja=1:no_images_in_sys
            % line_no+ja-no_images
             min_distotemp=sqrt((f(line_no+ja-no_images,1)-x_start+1-m)^2+(f(line_no+ja-no_images,2)-y_start+1-n)^2);
      
         if min_disto(line_no+ja-no_images)>min_distotemp
          
         min_disto(line_no+ja-no_images)=min_distotemp;
         closest_image_m(line_no+ja-no_images)=m+x_start;
         closest_image_n(line_no+ja-no_images)=n+y_start;
    
         end
         end
        end
         
    end
   end
end
end
clear imo s_check

line_no=line_no+1;

end %going through systems

if flag_lensingok==1
min_disto=min_disto*pix_scale;
rms=sqrt(sum((min_disto).^2)/length(f(:,1)))
final_chi2_pixels(ing)=sum((min_disto).^2/1^2)


else
    final_chi2_pixels(ing)=10^9
end
end

minchi=min(final_chi2_pixels);
indo_min=find(final_chi2_pixels==minchi);
display('best step was:')
indo_min
display('meaning a revised dls/ds factor of:')
k_system=(f(line_no-1,3)-minus+dg*(indo_min-1))
display('or a redshift of:')
calculate_z_for_k(z_l, z_s_main, k_system);
