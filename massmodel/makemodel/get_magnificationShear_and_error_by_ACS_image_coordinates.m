function [magn,magupvalue,magdownvalue,magupvalue_bypos,magdownvalue_bypos,gammanu,phinu,MJAn,MNAn]=get_magnificationShear_and_error_by_ACS_image_coordinates(xacs,yacs,zobject,zerr,zcluster,mainsystemz,location_err,image_to_plot_on)

%%put in the direct acs coordinates from the figur you have been working on
%%the location err should be taken usually as \pm0.5 arcsec in the position
%%of the critical curves meanning the value should be 10 pixels if acs
%%the script will out put both errors and you shall decide.
%%e.g.  [magn,magupvalue,magdownvalue,magupvalue_bypos,magdownvalue_bypos]=get_magnification_and_error_by_ACS_image_coordinates(2234,3170,7,1,0.33,2.63,10)
%%the script will show you the magnification image with a mark on, 
%%and when you click with the mouse on it it will show you the stiff image with a mark on 

[dlsratio]=calculate_k_for_system_out(zcluster, mainsystemz, zobject);
load transformed0201
%load transformed0201
alpha_x_ALL=alpha_x_ALL*dlsratio;
alpha_y_ALL=alpha_y_ALL*dlsratio;
%calculate magnification and mass:
for m=2:2998
   % m
 for  n=2:2998
     %
                da_x=(alpha_x_ALL(m+1,n)-alpha_x_ALL(m-1,n))/2.0;
                da_y=(alpha_y_ALL(m,n+1)-alpha_y_ALL(m,n-1))/2.0;
                da_x_dy=(alpha_x_ALL(m,n+1)-alpha_x_ALL(m,n-1))/2.0;
                da_y_dx=(alpha_y_ALL(m+1,n)-alpha_y_ALL(m-1,n))/2.0;
                poisson_ALL(m-1,n-1)=da_x+da_y;
                magnification_ALL(m-1,n-1)=abs(1.0/(1.0-poisson_ALL(m-1,n-1)+da_x*da_y-da_x_dy*da_y_dx));
                gammax(m-1,n-1)=0.5*(da_x-da_y);
                gammay(m-1,n-1)=da_y_dx;
                gamma_abs(m-1,n-1)=sqrt(gammax(m-1,n-1)^2+gammay(m-1,n-1)^2);
                phi_abs(m-1,n-1)=atan(gammay(m-1,n-1)/gammax(m-1,n-1));
                MJA(m-1,n-1)=1/(1-0.5*(poisson_ALL(m-1,n-1))-gamma_abs(m-1,n-1));
                MNA(m-1,n-1)=1/(1-0.5*(poisson_ALL(m-1,n-1))+gamma_abs(m-1,n-1));
                %magnification_with_sign(m,n)=(1.0/(1.0-poisson(m,n)+da_x*da_y-da_x_dy*da_y_dx));
                 
                end
end

magn=magnification_ALL(xacs-x_start,yacs-y_start)
magupvalue_bypos=max(max(magnification_ALL(xacs-x_start-location_err:xacs-x_start+location_err,yacs-y_start-location_err:yacs-y_start+location_err)));
magdownvalue_bypos=min(min(magnification_ALL(xacs-x_start-location_err:xacs-x_start+location_err,yacs-y_start-location_err:yacs-y_start+location_err)));
gammanu=gamma_abs(xacs-x_start,yacs-y_start)
%gammanupvalue_bypos=max(max(gamma_abs(xacs-x_start-location_err:xacs-x_start+location_err,yacs-y_start-location_err:yacs-y_start+location_err)));
%gammanvalue_bypos=min(min(gamma_abs(xacs-x_start-location_err:xacs-x_start+location_err,yacs-y_start-location_err:yacs-y_start+location_err)));
phinu=phi_abs(xacs-x_start,yacs-y_start)
MJAn=MJA(xacs-x_start,yacs-y_start)
MNAn=MNA(xacs-x_start,yacs-y_start)
MJAupvalue_bypos=max(max(MJA(xacs-x_start-location_err:xacs-x_start+location_err,yacs-y_start-location_err:yacs-y_start+location_err)))
MJAdownvalue_bypos=min(min(MJA(xacs-x_start-location_err:xacs-x_start+location_err,yacs-y_start-location_err:yacs-y_start+location_err)))
MNAupvalue_bypos=max(max(MNA(xacs-x_start-location_err:xacs-x_start+location_err,yacs-y_start-location_err:yacs-y_start+location_err)))
MNAdownvalue_bypos=min(min(MNA(xacs-x_start-location_err:xacs-x_start+location_err,yacs-y_start-location_err:yacs-y_start+location_err)))

[dlsratio]=calculate_k_for_system_out(zcluster, mainsystemz, zobject+zerr);
load transformed0201
alpha_x_ALL=alpha_x_ALL*dlsratio;
alpha_y_ALL=alpha_y_ALL*dlsratio;
%calculate magnification and mass:
for m=2:2998
   % m
 for  n=2:2998
     %
                da_x=(alpha_x_ALL(m+1,n)-alpha_x_ALL(m-1,n))/2.0;
                da_y=(alpha_y_ALL(m,n+1)-alpha_y_ALL(m,n-1))/2.0;
                da_x_dy=(alpha_x_ALL(m,n+1)-alpha_x_ALL(m,n-1))/2.0;
                da_y_dx=(alpha_y_ALL(m+1,n)-alpha_y_ALL(m-1,n))/2.0;
                poisson_ALL(m-1,n-1)=da_x+da_y;
                magnification_ALL(m-1,n-1)=abs(1.0/(1.0-poisson_ALL(m-1,n-1)+da_x*da_y-da_x_dy*da_y_dx));
                %magnification_with_sign(m,n)=(1.0/(1.0-poisson(m,n)+da_x*da_y-da_x_dy*da_y_dx));
                 gammax(m-1,n-1)=0.5*(da_x-da_y);
                gammay(m-1,n-1)=da_y_dx;
                gamma_abs(m-1,n-1)=sqrt(gammax(m-1,n-1)^2+gammay(m-1,n-1)^2);
                phi_abs(m-1,n-1)=atan(gammay(m-1,n-1)/gammax(m-1,n-1));
              MJA(m-1,n-1)=1/(1-0.5*(poisson_ALL(m-1,n-1))-gamma_abs(m-1,n-1));
                MNA(m-1,n-1)=1/(1-0.5*(poisson_ALL(m-1,n-1))+gamma_abs(m-1,n-1));
               
                end
end

magupvalue=magnification_ALL(xacs-x_start,yacs-y_start)
gammaupvalue=gamma_abs(xacs-x_start,yacs-y_start)
phiupvalue=phi_abs(xacs-x_start,yacs-y_start)
MJAnupvalue=MJA(xacs-x_start,yacs-y_start)
MNAnupvalue=MNA(xacs-x_start,yacs-y_start)
[dlsratio]=calculate_k_for_system_out(zcluster, mainsystemz, zobject-zerr);
load transformed0201
alpha_x_ALL=alpha_x_ALL*dlsratio;
alpha_y_ALL=alpha_y_ALL*dlsratio;
%calculate magnification and mass:
for m=2:2998
   % m
 for  n=2:2998
     %
                da_x=(alpha_x_ALL(m+1,n)-alpha_x_ALL(m-1,n))/2.0;
                da_y=(alpha_y_ALL(m,n+1)-alpha_y_ALL(m,n-1))/2.0;
                da_x_dy=(alpha_x_ALL(m,n+1)-alpha_x_ALL(m,n-1))/2.0;
                da_y_dx=(alpha_y_ALL(m+1,n)-alpha_y_ALL(m-1,n))/2.0;
                poisson_ALL(m-1,n-1)=da_x+da_y;
                magnification_ALL(m-1,n-1)=abs(1.0/(1.0-poisson_ALL(m-1,n-1)+da_x*da_y-da_x_dy*da_y_dx));
                %magnification_with_sign(m,n)=(1.0/(1.0-poisson(m,n)+da_x*da_y-da_x_dy*da_y_dx));
                gammax(m-1,n-1)=0.5*(da_x-da_y);
                gammay(m-1,n-1)=da_y_dx;
                gamma_abs(m-1,n-1)=sqrt(gammax(m-1,n-1)^2+gammay(m-1,n-1)^2);
                phi_abs(m-1,n-1)=atan(gammay(m-1,n-1)/gammax(m-1,n-1));
               MJA(m-1,n-1)=1/(1-0.5*(poisson_ALL(m-1,n-1))-gamma_abs(m-1,n-1));
                MNA(m-1,n-1)=1/(1-0.5*(poisson_ALL(m-1,n-1))+gamma_abs(m-1,n-1));
               
                end
end
magdownvalue=magnification_ALL(xacs-x_start,yacs-y_start)
gammadonvalue=gamma_abs(xacs-x_start,yacs-y_start)
phidownvalue=phi_abs(xacs-x_start,yacs-y_start)
magnification_ALL(xacs-x_start-50:xacs-x_start+50,yacs-y_start-50:yacs-y_start+50)=3000;
MJAndownvalue=MJA(xacs-x_start,yacs-y_start)
MNAndownvalue=MNA(xacs-x_start,yacs-y_start)

magnification_im_ALL=magnification_ALL'; 
imagesc(log(magnification_im_ALL+1)); set(gca,'YDir','normal');
display('press with mouse on the image to continue');
waitforbuttonpress
toring=log((magnification_im_ALL+1));
fl=imread(image_to_plot_on);
for m=1:(x_size-i_x-5)
   % m
 for  n=1:(y_size-i_y-5)
     
     if (toring(m,n)>=6)
         
         fl(3000-m,n,:)=255;
     end

 end
end

%fl=imread(image_to_plot_on);
% fl((3000-(yacs-y_start-10)):(3000-(yacs-y_start+10)),xacs-x_start-10:xacs-x_start+10,1)=255;
% fl((3000-(yacs-y_start-10)):(3000-(yacs-y_start+10)),xacs-x_start-10:xacs-x_start+10,2)=0;
% fl((3000-(yacs-y_start-10)):(3000-(yacs-y_start+10)),xacs-x_start-10:xacs-x_start+10,3)=0;
imshow(fl);
%%plot location on stiff image