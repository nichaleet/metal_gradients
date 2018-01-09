im=fitsread('cswa165_SDSSJ0105+0144_deflections.fits');
srcx=im(:,:,2);
srcy=im(:,:,3)  ;
[leng_y,leng_x]=size(srcx);
[imx,imy] = meshgrid(1:leng_x,1:leng_y);
a_x = srcx-imx;
a_y = srcy-imy;
[da_x_dy,da_x_dx]=gradient(flipud(a_x));
[da_y_dy,da_y_dx]=gradient(flipud(a_y));
poisson_ALL=da_x_dx+da_y_dy;
magnification_ALL=abs(1./(1-poisson_ALL+da_x_dx.*da_y_dy-da_x_dy.*da_y_dx));
fitswrite(flipud(magnification_ALL),'cswa165_magnification.fits');
imref=fitsread('/scr2/nichal/workspace/output/cswa165_Ha_tlc_Kn2_handmosaic_sky_2hr_acube.fits');
detection=imref(:,:,5);
detection=flipud(detection);
detection=detection(1:75,1:66);
%detection = detection(:,40:66);
%magnification_ALL = magnification_ALL(:,40:66);
figure(1)
image(magnification_ALL*10.);
figure(2)
image(detection);
good=find(isfinite(detection) & detection>10 & magnification_ALL>1.);
mag=mean(magnification_ALL(good))

