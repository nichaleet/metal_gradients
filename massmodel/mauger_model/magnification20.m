im=fitsread('cswa20_SDSSJ1441+1441_deflections.fits');
imx=im(:,:,2);
imy=im(:,:,3);
[da_x_dy,da_x_dx]=gradient(flipud(imx));
[da_y_dy,da_y_dx]=gradient(flipud(imy));
poisson_ALL=da_x_dx+da_y_dy;
magnification_ALL=abs(1./(1-poisson_ALL+da_x_dx.*da_y_dy-da_x_dy.*da_y_dx));
fitswrite(flipud(magnification_ALL),'cswa20_magnification.fits');
imref=fitsread('/scr2/nichal/workspace/output/cswa20_Ha_tlc_Hn2_handmosaic_scaledsky_1hr_acube.fits');
detection=imref(:,:,5);
detection=flipud(detection);
detection=detection(1:71,1:66)
detection = detection(:,40:66);
magnification_ALL = magnification_ALL(:,40:66);
figure(1)
image((magnification_ALL*10))
figure(2)
image(detection)
good=find(isfinite(detection) & detection>10);
mag=mean(magnification_ALL(good))

