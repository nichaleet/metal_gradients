function [final_chi2_pixels]=go_try_specific_mc_gaussianBCGSplinerms(q,s,k_new,k_gal,gamma,phi,BCG,galaxies_file,images_file,pix_scale)
%%here I don;t save the files!
%[alpha_x,alpha_y] = cluster_lense_model_loops('model_colors_1206new.txt',171,q);
file_data=load(galaxies_file);
%-------------------------------------
%note: BCG must be first line!!!!
%------------------------------------

%i_x and i_y represent the starting x and y in pixels of the desired frame
%given in the coordinate system of the galaxies coordinate file: 
% min_x=min(file_data(:,2));
% max_x=max(file_data(:,2));
% min_y=min(file_data(:,3));
% max_y=max(file_data(:,3));
i_x=1300;
i_y=1300;
leng=3000;
x_start=i_x;
y_start=i_y;
m=i_x;
n=i_y;
obj=1;

x_size=4300;%=3916-500=3416
y_size=4300;%=4279-500=3779
%initializing the alpha matrices:
 alpha_x_gpu(1:x_size-i_x+1,1:y_size-i_y+1)=0;
 alpha_y_gpu(1:x_size-i_x+1,1:y_size-i_y+1)=0;
 %display('cgpuArrayTransfer')
%  alpha_x_gpu=gpuArray(alpha_x_ram);
%  alpha_y_gpu=gpuArray(alpha_y_ram);
 clear alpha_x_ram alpha_y_ram
 display('calculating galaxies field')
[YI,XI]=meshgrid(1:y_size-i_y+1,1:x_size-i_x+1);
% XI=gpuArray(XI_ram);
% YI=gpuArray(YI_ram);
clear XI_ram YI_ram
XI=XI+i_x-1;
YI=YI+i_y-1;
obj=1;
    dif_x=-(file_data(obj,2)-XI);
    dif_y=-(file_data(obj,3)-YI);
    shoresh=dif_x.^2+dif_y.^2;
    if file_data(obj,2)==round(file_data(obj,2)) && file_data(obj,3)==round(file_data(obj,3))
        shoresh(file_data(obj,2)-i_x+1,file_data(obj,3)-i_y+1)=1.0;
    end
    alpha_x_gpu=alpha_x_gpu+BCG*file_data(obj,7)*dif_x./shoresh.^(0.5*q);
    alpha_y_gpu=alpha_y_gpu+BCG*file_data(obj,7)*dif_y./shoresh.^(0.5*q);
    
for obj=2:length(file_data(:,1))
    
    dif_x=-(file_data(obj,2)-XI);
    dif_y=-(file_data(obj,3)-YI);
    shoresh=dif_x.^2+dif_y.^2;
    if file_data(obj,2)==round(file_data(obj,2)) && file_data(obj,3)==round(file_data(obj,3))
        shoresh(file_data(obj,2)-i_x+1,file_data(obj,3)-i_y+1)=1.0;
    end
    alpha_x_gpu=alpha_x_gpu+file_data(obj,7)*dif_x./shoresh.^(0.5*q);
    alpha_y_gpu=alpha_y_gpu+file_data(obj,7)*dif_y./shoresh.^(0.5*q);
end
%[magnification]=calculate_magnification;
factor=1.8e-1;
%calculate magnification and mass:
alpha_x=alpha_x_gpu;
alpha_y=alpha_y_gpu;
%display('starting magnif first')
[da_x_dy,da_x_dx]=gradient(factor*alpha_x);
[da_y_dy,da_y_dx]=gradient(factor*alpha_y);
poisson=da_x_dx+da_y_dy;
%smooth_comp=poisson;
%magnification=abs(1./(1-poisson+da_x_dx.*da_y_dy-da_x_dy.*da_y_dx));
%display('ended magnif first')
%[smooth_comp]=smooth_DM(s);

% h=fspecial('gaussian',2500,s);
% smooth_comp=imfilter(poisson,h);
X=(1:(x_size-i_x-1)/5+1);
Y=(1:(y_size-i_y-1)/5+1);
% size(X)
% size(Y)
c(1:(x_size-i_x-1)/5+1,1:(y_size-i_y-1)/5+1)=1;
size(c)
size(poisson)
for i=1:(x_size-i_x-1)/5+1
    for j=1:(y_size-i_y-1)/5+1
        poisson_new(i,j)=poisson(i*5-1,j*5-1);
    end
end

%clear poisson
% size(poisson_new)

%set degree of polynom
ZI=polyfitweighted2(Y,X,poisson_new,s,c);
smooth_comp_first=polyval2(ZI,Y,X);
[XI, YI]=meshgrid(1:0.2:(x_size-i_x-1)/5+1);
[m,n]=size(poisson_new);
smooth_comp=interp2(1:m,1:n,smooth_comp_first,XI,YI,'*spline');


%[alpha_x_DM,alpha_y_DM] = cluster_lense_model_loops_DM;

display('started fourier DM')
x_factor(1:2*length(smooth_comp(:,1)),1:2*length(smooth_comp(1,:)))=0;
y_factor(1:2*length(smooth_comp(:,1)),1:2*length(smooth_comp(1,:)))=0;
%g_factor(1:2*length(smooth_comp(:,1)),1:2*length(smooth_comp(1,:)))=0;
[Y,X]=meshgrid(-length(smooth_comp(:,1)):length(smooth_comp(:,1))-1,-length(smooth_comp(1,:)):length(smooth_comp(1,:))-1);
x_factor=X./(X.^2+Y.^2);
y_factor=Y./(X.^2+Y.^2);
%g_factor=exp(-(X.^2+Y.^2)./(2.0*s*s))./(2.0*pi*s*s);
x_factor_shift=fftshift(x_factor);
y_factor_shift=fftshift(y_factor);
%g_factor_shift=fftshift(g_factor);
x_factor_shift(1,1)=0;
y_factor_shift(1,1)=0;
xfourier=fft2(x_factor_shift);
yfourier=fft2(y_factor_shift);
%gfourier=fft2(g_factor_shift);

new_smooth(1:2*length(smooth_comp(:,1)), 1:2*length(smooth_comp(1,:)))=0;%mean2(smooth_comp_first);
new_smooth(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)))=smooth_comp;

fourier_smooth_new=fft2(new_smooth);

alpha_x_DMf=ifft2(xfourier.*fourier_smooth_new);
alpha_y_DMf=ifft2(yfourier.*fourier_smooth_new);
%clear new_smooth xfourier yfourier
alpha_x_DMfs=alpha_x_DMf(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)));
alpha_y_DMfs=alpha_y_DMf(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)));
clear alpha_x_DMf alpha_y_DMf
% [XI, YI]=meshgrid(1:0.2:(x_size-i_x-1)/5+1);
% [m,n]=size(poisson_new);
% alpha_x_DM=interp2(1:m,1:n,alpha_y_DMfs,XI,YI,'*spline');
% alpha_y_DM=interp2(1:m,1:n,alpha_x_DMfs,XI,YI,'*spline');
alpha_x_DM=alpha_x_DMfs/2520;
alpha_y_DM=alpha_y_DMfs/2520;
%save 'alpha_xtry.txt' alpha_x_DM -ascii
clear  alpha_x_DMfs alpha_y_DMfs
display('finished DM fourier')

%[magnification_ALL]=calculate_magnification_ALL_general(k_new,k_gal,gamma,phi);

rat=900/alpha_x_DM(x_size-i_x-5,round((x_size-i_x)/2));
rat_gal=10000/alpha_x(x_size-i_x-5,round((x_size-i_x)/2));

factor=1.8*k_new;

gammacos2phi=gamma*cos(2*phi);
gammasin2phi=gamma*sin(2*phi);
[YI,XI]=meshgrid(1:(y_size-i_y),1:(x_size-i_x));
alpha_x_add=gammacos2phi*XI+gammasin2phi*YI;
alpha_y_add=gammasin2phi*XI-gammacos2phi*YI;

%save tillnow1
alpha_x_ALL=(k_gal*alpha_x(1:leng,1:leng)*rat_gal+(1-k_gal)*alpha_x_DM(1:leng,1:leng)*rat)*factor+alpha_x_add(1:leng,1:leng);%+k_gal*k_new*k_plus*alpha_x_plus(1:leng,1:leng);
alpha_y_ALL=(k_gal*alpha_y(1:leng,1:leng)*rat_gal+(1-k_gal)*alpha_y_DM(1:leng,1:leng)*rat)*factor+alpha_y_add(1:leng,1:leng);%+k_gal*k_new*k_plus*alpha_y_plus(1:leng,1:leng);
%this factor should be itterated. Always with plus sign: 
display('finished calculating total def field')
%calculate magnification and mass:
[da_x_dy,da_x_dx]=gradient(alpha_x_ALL);
[da_y_dy,da_y_dx]=gradient(alpha_y_ALL);
poisson_ALL=da_x_dx+da_y_dy;
magnification_ALL=abs(1./(1-poisson_ALL+da_x_dx.*da_y_dy-da_x_dy.*da_y_dx));
%magnification_ALL_with_sign=(1./(1-poisson_ALL+da_x_dx.*da_y_dy-da_x_dy.*da_y_dx));
clear alpha_x_DM alpha_x_add alpha_y_DM alpha_y_add X Y XI YI da_x_dx da_y_dy da_x_dy da_y_dx dif_x dif_y poisson
clear x_factor y_factor x_factor_shift y_factor_shift
save bestmodel_aCovtry_newallsysV1
%display('finished calculating tmagnification and kappa')
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

    imo(a-3:a+3,b-3:b+3)=100;
 
for m=a-3:a+3
 for n=b-3:b+3
     im_x(countr)=m-round(alpha_x_ALL(m,n)*f(line_no,3));
     im_y(countr)=n-round(alpha_y_ALL(m,n)*f(line_no,3)) ;
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
s_check(round(imxavg)-3:round(imxavg)+3,round(imyavg)-3:round(imyavg)+3)=100;
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
   
    im_x=m-round(alpha_x_ALL(m,n)*f(line_no,3));
    im_y=n-round(alpha_y_ALL(m,n)*f(line_no,3));
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
final_chi2_pixels=sum((min_disto).^2/1.9^2)


% final_chi2_pixels=sum(sum_disto);
% final_chi2_average_of_squares=sum(sum_disto_square);
% %now calculate all quantities
% avg_per_all_images=sum(ave_per_image(:))/no_systems * pix_scale
% display('arcseconds')
% rms=sqrt((sum(ave_per_image(:).^2))/no_systems) * pix_scale
% display('arcseconds')
% display('this is chi2 by square of average:')
% final_chi2_pixels
% display('this is chi2 by average of square:')
% final_chi2_average_of_squares
% display('reduced chi2 by square of average:')
% final_chi2_pixels/abs(length(f(:,1))-2)
% display('reduced chi2 by average of square:')
% final_chi2_average_of_squares/abs(length(f(:,1))-2)
else
    final_chi2_pixels=10^9
end
