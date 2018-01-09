function MCstep_PowerLaw_and_GaussianOrSplineByEye


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for Zitrin & Broadhurst lens modeling script version 06.12.2012
%cite (Zitrin et al., 2009, MNRAS, 396, 1985) 
%Written by Adi Zitrin with Tom Broadhurst 
%Suggested and implemented improvements by other contributors is highly
%appreciated especially Matthias Bartelmann, Gregor Seidel, Irene Sendra

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function is creating a mass model and calculating the chi2 and rms 
%in each step of the MCMC
%example:
%[alpha_x,alpha_y] = cluster_lense_model_loops('model_colors_1206new.txt',171,q);
%(file_data is model galaxies file)

%NOTE: user has to edit one small section below

%define manually:
load MCinput
params_data=[1.3,20,0.13,0.15,0,0,0,0,0,0,0,0,1,1,1];




 %load FinalChain
  %   params_data=best_vec(1:Nparams);
  % params_data
%good old with 1.3 fixed and 20 free gals rms~6.2 chi2~16.35: [1.30000000000000,10,0.596246934026259,0.0196268879175615,0.113296641670878,0.565162798633010,1.92785703361407,0.140901566980991,1.39947118906233,1.36980409920663,0.385116014525295,0.420126183356900,0.232554247353881,0.0116932369044308,0.0566790064916258,0.841237815142694,0.167286011641672,0.0797677864052256,0.266965660533915,0.134229535291840,0.197420256969290,0.0253172155893736,0.398451419408569,0.474430761549718,0.0656104041779518,0.214333754707173;]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%persistent fwl resto file_data leng2 x_factor y_factor x_factor_shift y_factor_shift xfourier yfourier    x_start y_start x_size y_size X1 Y1 XI2 YI2 grid_image
%persistent Xg Yg XI2gpu YI2gpu alpha_xf_ram YI_ram XI_ram  YI_ram2 XI_ram2 XI YI X Y f 
%persistent arcsecscale arcsecscalecm lumidist1 ignore1 ignore2 lumidist
%if isempty(file_data)
%    arcsecscale= tan(1/3600*pi/180)*ad_dist(cluster_redshift);% pc per arcsec;
%arcsecscalecm=arcsecscale* 3.08568025*10^18 ; %cm per arcsec
%[lumidist1,ignore1,ignore2]=lum_dist(cluster_redshift);
%lumidist=lumidist1/ arcsecscale*pix_scale; %cm per arcsec
    resto=1%Wresolution;
        fwl=shapes_file; 
    %jumpi=jumpi*10/Wresolution;
file_data=galaxies_file;
x_start=i_x;
y_start=i_y;
x_size=i_x+leng;%=3916-500=3416
y_size=i_y+leng;%=4279-500=3779
leng2=round(leng/resto);

%initializing the alpha matrices:
alpha_xf_ram(1:3*leng2, 1:3*leng2)=0;
[YI_ram,XI_ram]=meshgrid(1:resto:y_size-i_y+1,1:resto:x_size-i_x+1);
if strcmp(GPU_OR_CPU,'GPU')==1
XI=gpuArray(XI_ram);
YI=gpuArray(YI_ram);
elseif strcmp(GPU_OR_CPU,'CPU')==1
XI=XI_ram;
YI=YI_ram;  
clear XI_ram YI_ram 
end

XI=XI+i_x-1;
YI=YI+i_y-1;
[Y,X]=meshgrid(-leng:resto:leng-1,-leng:resto:leng-1);
[Yg,Xg]=meshgrid(-1.5*leng:resto:1.5*leng-1,-1.5*leng:resto:1.5*leng-1);
X1=(1:(x_size-i_x-1)/(resto)+1);
Y1=(1:(y_size-i_y-1)/(resto)+1);
[YI2gpu,XI2gpu]=meshgrid(1:resto:(y_size-i_y),1:resto:(x_size-i_x));
if strcmp(GPU_OR_CPU,'GPU')==1
XI2=gpuArray(XI2gpu);
YI2=gpuArray(YI2gpu);
elseif strcmp(GPU_OR_CPU,'CPU')==1
XI2=XI2gpu;
YI2=YI2gpu;
clear XI2gpu YI2gpu
end

f=images_file;
[YI_ram2,XI_ram2]=meshgrid(1:y_size-i_y,1:x_size-i_x);
if strcmp(GPU_OR_CPU,'GPU')==1
grid_image = gpuArray(complex(XI_ram2,YI_ram2));
elseif strcmp(GPU_OR_CPU,'CPU')==1
grid_image = (complex(XI_ram2,YI_ram2));
end
grid_image=grid_image(1:resto:leng,1:resto:leng);
x_factor=X./(X.^2+Y.^2);
y_factor=Y./(X.^2+Y.^2);
x_factor_shift=fftshift(x_factor);
y_factor_shift=fftshift(y_factor);
x_factor_shift(1,1)=0;
y_factor_shift(1,1)=0;
xfourier=fft2(x_factor_shift);
yfourier=fft2(y_factor_shift);

%if strcmp(SplineOrGaussian,'Gaussian')==1

%end

%end
%resto=Wresolution;
%-------------------------------------
%note: q,s,k_new,k_gal,gamma,phi,BCG should be the order of vector of input
%params
%------------------------------------
%='Gaussian';
ModeRun=1;
%renormalization of the def field according to the power law
%preliminary expression should be revised! :
% if params_data(1)>1
% params_data(4)=params_data(4)/3 * ((params_data(1)-1)*10+params_data(1)^2)^params_data(1);
% elseif params_data(1)==1
% params_data(4)=params_data(4)/3 ;
%end
jumpi=60*10/resto;
q=params_data(1);
s=params_data(2);
k_new=params_data(3);
k_gal=params_data(4);
gamma=params_data(5);
phi=params_data(6);
rci=params_data(7);
BCGel=params_data(8);
BCGPA=params_data(9);
rci2=params_data(10);
BCGel2=params_data(11);
BCGPA2=params_data(12);
rci
BCGel
BCGPA
rci2
BCGel2
BCGPA2
%params_data(1:6)
% if ModeRun==3
if CosmoFree==1
    omegaM=params_data(13);
    w_0=params_data(14);
    omegaL=params_data(15);               %flat
    wa=params_data(16);
    H_zero=params_data(17);
end
% if ModeRun==3
% wa=params_data(6+No_of_Bright_Free+No_of_redshift_Free+5);
% %wa
% w_0=params_data(6+No_of_Bright_Free+No_of_redshift_Free+4);
% %w_0
% omegaL=params_data(6+No_of_Bright_Free+No_of_redshift_Free+3);
% %omegaL
% omegaM=params_data(6+No_of_Bright_Free+No_of_redshift_Free+2);
% %omegaM
% H_zero=params_data(6+No_of_Bright_Free+No_of_redshift_Free+1);
% %H_zero
% end
%-------------------------------------------------------------------------------

%mean2(smooth_comp_first);
%file_data1=galaxies_file;
% if No_of_Bright_Free>=1
% for i=1:No_of_Bright_Free
% maxtemp=max(file_data1(:,7));%7th column is the flux
% Ind_Free_Gal(i)=find(file_data1(:,7)==maxtemp);
% file_data1(Ind_Free_Gal(i),:)=0;
% end
% elseif No_of_Bright_Free==0
%    Ind_Free_Gal(1)=0;
% end
% 
% %Ind_Free_Gal(1:10)=[187 83 142 159 16 139 189 118 36 110];
% Ind_Free_Gal

if strcmp(GPU_OR_CPU,'GPU')==1
alpha_xf=gpuArray(alpha_xf_ram);
elseif strcmp(GPU_OR_CPU,'CPU')==1
    alpha_xf=(alpha_xf_ram);
end
 t=1;
        epsgal=BCGel;%0.1765;
         PAgal=BCGPA*pi/180;%pi/2+42.5*pi/180;
         %  epsgal=0.228;
         %PAgal=pi/2+45*pi/180;
    %xcent=(file_data(t,2)-i_x)+1;%file_data(t,2)-i_x+1;
      % ycent=(file_data(t,3)-i_y)+1;%file_data(t,3)-i_y+1;
        bcg_x=(file_data(t,2)-i_x)+1;%file_data(t,2)-i_x+1;
      bcg_y=(file_data(t,3)-i_y)+1;%file_data(t,3)-i_y+1;
       delta_x1=(Xg)*cos(PAgal)+(Yg)*sin(PAgal);%*cos(phi);
       delta_y1=-(Xg)*sin(PAgal)+(Yg)*cos(PAgal);%*sin(phi);
        
      rtemp=sqrt(1/(1+epsgal)^2*delta_x1.^2+1/(1-epsgal)^2*delta_y1.^2);
        
    if file_data(t,2)>i_x+1-leng2 && file_data(t,2)<i_x+leng-1+leng2 && file_data(t,3)>i_y+1-leng2 && file_data(t,3)<i_y+leng-1+leng2 
    if find(Ind_Free_Gal(:)==t)>0
        atempto=find(Ind_Free_Gal(:)==t);
        if CosmoFree==1
        alpha_xf(round((file_data(t,2)-i_x+leng)/resto),round((file_data(t,3)-i_y+leng)/resto))=file_data(t,7)*params_data(atempto+17);%%need to add weights
        else
            alpha_xf(round((file_data(t,2)-i_x+leng)/resto),round((file_data(t,3)-i_y+leng)/resto))=file_data(t,7)*params_data(atempto+12);%%need to add weights
        end
    else
        alpha_xf(round((file_data(t,2)-i_x+leng)/resto),round((file_data(t,3)-i_y+leng)/resto))=file_data(t,7);%%need to add weights
    
    end
    end

alpha_yf=alpha_xf;

xf_factor=Xg./((rtemp).^2 +rci^2).^(0.5*q);
yf_factor=Yg./((rtemp).^2 +rci^2).^(0.5*q);

xf_factor_shift=fftshift(xf_factor);
yf_factor_shift=fftshift(yf_factor);

xf_factor_shift(1,1)=0;
yf_factor_shift(1,1)=0;
xffourier=fft2(xf_factor_shift);
yffourier=fft2(yf_factor_shift);
fourier_alpha_xf=fft2(alpha_xf);
fourier_alpha_yf=fft2(alpha_yf);

alpha_x_ff=ifft2(xffourier.*fourier_alpha_xf);
alpha_y_ff=ifft2(yffourier.*fourier_alpha_yf);
%alpha_y_DMf=ifft2(yfourier.*fourier_smooth_new);
%clear new_smooth xfourier yfourier
alpha_x_ffs=alpha_x_ff(leng2+1:leng2+leng2,leng2+1:leng2+leng2);
alpha_y_ffs=alpha_y_ff(leng2+1:leng2+leng2,leng2+1:leng2+leng2);
%alpha_y_DMfs=alpha_y_DMf(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)));
clear alpha_x_ff alpha_y_ff
if strcmp(GPU_OR_CPU,'GPU')==1
alpha_x1=real(gather(alpha_x_ffs));%./lumidist^(1-params_data(1)).*(2-params_data(1))/params_data(1)^2.2;
alpha_y1=real(gather(alpha_y_ffs));%./lumidist^(1-params_data(1)).*(2-params_data(1))/params_data(1)^2.2;
elseif strcmp(GPU_OR_CPU,'CPU')==1
alpha_x1=real((alpha_x_ffs));%./lumidist^(1-params_data(1)).*(2-params_data(1))/params_data(1)^2.2;
alpha_y1=real((alpha_y_ffs));%./lumidist^(1-params_data(1)).*(2-params_data(1))/params_data(1)^2.2;
end
% 
 clear alpha_xf

if strcmp(GPU_OR_CPU,'GPU')==1
alpha_xf=gpuArray(alpha_xf_ram);
elseif strcmp(GPU_OR_CPU,'CPU')==1
    alpha_xf=(alpha_xf_ram);
end
 t=2;
        epsgal=BCGel2;%0.1765;
         PAgal=BCGPA2*pi/180;%pi/2+42.5*pi/180;
         %  epsgal=0.228;
         %PAgal=pi/2+45*pi/180;
    %xcent=(file_data(t,2)-i_x)+1;%file_data(t,2)-i_x+1;
      % ycent=(file_data(t,3)-i_y)+1;%file_data(t,3)-i_y+1;
      
       delta_x1=(Xg)*cos(PAgal)+(Yg)*sin(PAgal);%*cos(phi);
       delta_y1=-(Xg)*sin(PAgal)+(Yg)*cos(PAgal);%*sin(phi);
        
      rtemp=sqrt(1/(1+epsgal)^2*delta_x1.^2+1/(1-epsgal)^2*delta_y1.^2);
        
    if file_data(t,2)>i_x+1-leng2 && file_data(t,2)<i_x+leng-1+leng2 && file_data(t,3)>i_y+1-leng2 && file_data(t,3)<i_y+leng-1+leng2 
    if find(Ind_Free_Gal(:)==t)>0
        atempto=find(Ind_Free_Gal(:)==t);
        if CosmoFree==1
        alpha_xf(round((file_data(t,2)-i_x+leng)/resto),round((file_data(t,3)-i_y+leng)/resto))=file_data(t,7)*params_data(atempto+17);%%need to add weights
        else
            alpha_xf(round((file_data(t,2)-i_x+leng)/resto),round((file_data(t,3)-i_y+leng)/resto))=file_data(t,7)*params_data(atempto+12);%%need to add weights
        end
    else
        alpha_xf(round((file_data(t,2)-i_x+leng)/resto),round((file_data(t,3)-i_y+leng)/resto))=file_data(t,7);%%need to add weights
    
    end
    end

alpha_yf=alpha_xf;

xf_factor=Xg./((rtemp).^2 +rci2^2).^(0.5*q);
yf_factor=Yg./((rtemp).^2 +rci2^2).^(0.5*q);

xf_factor_shift=fftshift(xf_factor);
yf_factor_shift=fftshift(yf_factor);

xf_factor_shift(1,1)=0;
yf_factor_shift(1,1)=0;
xffourier=fft2(xf_factor_shift);
yffourier=fft2(yf_factor_shift);
fourier_alpha_xf=fft2(alpha_xf);
fourier_alpha_yf=fft2(alpha_yf);

alpha_x_ff=ifft2(xffourier.*fourier_alpha_xf);
alpha_y_ff=ifft2(yffourier.*fourier_alpha_yf);
%alpha_y_DMf=ifft2(yfourier.*fourier_smooth_new);
%clear new_smooth xfourier yfourier
alpha_x_ffs=alpha_x_ff(leng2+1:leng2+leng2,leng2+1:leng2+leng2);
alpha_y_ffs=alpha_y_ff(leng2+1:leng2+leng2,leng2+1:leng2+leng2);
%alpha_y_DMfs=alpha_y_DMf(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)));
clear alpha_x_ff alpha_y_ff
if strcmp(GPU_OR_CPU,'GPU')==1
alpha_x2=real(gather(alpha_x_ffs));%./lumidist^(1-params_data(1)).*(2-params_data(1))/params_data(1)^2.2;
alpha_y2=real(gather(alpha_y_ffs));%./lumidist^(1-params_data(1)).*(2-params_data(1))/params_data(1)^2.2;
elseif strcmp(GPU_OR_CPU,'CPU')==1
alpha_x2=real((alpha_x_ffs));%./lumidist^(1-params_data(1)).*(2-params_data(1))/params_data(1)^2.2;
alpha_y2=real((alpha_y_ffs));%./lumidist^(1-params_data(1)).*(2-params_data(1))/params_data(1)^2.2;
end
% 
 clear alpha_xf 
 
if strcmp(GPU_OR_CPU,'GPU')==1
alpha_xf=gpuArray(alpha_xf_ram);
elseif strcmp(GPU_OR_CPU,'CPU')==1
    alpha_xf=(alpha_xf_ram);
end

for t=3:length(file_data(:,1))
    if file_data(t,2)>i_x+1-leng2 && file_data(t,2)<i_x+leng-1+leng2 && file_data(t,3)>i_y+1-leng2 && file_data(t,3)<i_y+leng-1+leng2 
    if find(Ind_Free_Gal(:)==t)>0
        atempto=find(Ind_Free_Gal(:)==t);
       if CosmoFree==1
        alpha_xf(round((file_data(t,2)-i_x+leng)/resto),round((file_data(t,3)-i_y+leng)/resto))=file_data(t,7)*params_data(atempto+17);%%need to add weights
        else
            alpha_xf(round((file_data(t,2)-i_x+leng)/resto),round((file_data(t,3)-i_y+leng)/resto))=file_data(t,7)*params_data(atempto+12);%%need to add weights
       end
    else
        alpha_xf(round((file_data(t,2)-i_x+leng)/resto),round((file_data(t,3)-i_y+leng)/resto))=file_data(t,7);%%need to add weights
    
    end
    end
end
alpha_yf=alpha_xf;

xf_factor=Xg./(Xg.^2+Yg.^2).^(0.5*q);
yf_factor=Yg./(Xg.^2+Yg.^2).^(0.5*q);

xf_factor_shift=fftshift(xf_factor);
yf_factor_shift=fftshift(yf_factor);

xf_factor_shift(1,1)=0;
yf_factor_shift(1,1)=0;
xffourier=fft2(xf_factor_shift);
yffourier=fft2(yf_factor_shift);
fourier_alpha_xf=fft2(alpha_xf);
fourier_alpha_yf=fft2(alpha_yf);

alpha_x_ff=ifft2(xffourier.*fourier_alpha_xf);
alpha_y_ff=ifft2(yffourier.*fourier_alpha_yf);
%alpha_y_DMf=ifft2(yfourier.*fourier_smooth_new);
%clear new_smooth xfourier yfourier
alpha_x_ffs=alpha_x_ff(leng2+1:leng2+leng2,leng2+1:leng2+leng2);
alpha_y_ffs=alpha_y_ff(leng2+1:leng2+leng2,leng2+1:leng2+leng2);
%alpha_y_DMfs=alpha_y_DMf(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)));
clear alpha_x_ff alpha_y_ff
if strcmp(GPU_OR_CPU,'GPU')==1
alpha_x=real(gather(alpha_x_ffs));%./lumidist^(1-params_data(1)).*(2-params_data(1))/params_data(1)^2.2;
alpha_y=real(gather(alpha_y_ffs));%./lumidist^(1-params_data(1)).*(2-params_data(1))/params_data(1)^2.2;
elseif strcmp(GPU_OR_CPU,'CPU')==1
alpha_x=real((alpha_x_ffs));%./lumidist^(1-params_data(1)).*(2-params_data(1))/params_data(1)^2.2;
alpha_y=real((alpha_y_ffs));%./lumidist^(1-params_data(1)).*(2-params_data(1))/params_data(1)^2.2;
end
alpha_x=alpha_x+alpha_x1+alpha_x2;
alpha_y=alpha_y+alpha_y1+alpha_y2;
% 
%image(alpha_x_ffs)
% alpha_x_DM=alpha_x_DMfs/2520;
% alpha_y_DM=alpha_y_DMfs/2520;
% 
% clear  alpha_x_DMfs alpha_y_DMfs
%------------------------------------------------------------------------------

%[magnification]=calculate_magnification;

%calculate magnification and mass:
% alpha_x=gather(alpha_x_gpu);
% alpha_y=gather(alpha_y_gpu);
%display('starting magnif first')
[da_x_dy,da_x_dx]=gradient(alpha_x/resto);
[da_y_dy,da_y_dx]=gradient(alpha_y/resto);
poisson=(da_x_dx+da_y_dy)/2;
image(fliplr(poisson*2/10000)')
waitforbuttonpress
%smooth_comp=poisson;
%magnification=abs(1./(1-poisson+da_x_dx.*da_y_dy-da_x_dy.*da_y_dx));
%display('ended magnif first')
%[smooth_comp]=smooth_DM(s);

% h=fspecial('gaussian',2500,s);
% smooth_comp=imfilter(poisson,h);
if strcmp(SplineOrGaussian,'Spline')==1


% size(X)
% size(Y)


% for i=1:(x_size-i_x-1)/(2*2)+1
%     for j=1:(y_size-i_y-1)/(2*2)+1
      %  poisson((1:round((x_size-i_x-1)/(2))),(1:round((y_size-i_y-1)/(2))))=poisson( (1:1:round((x_size-i_x-1))),(1:1:round((y_size-i_y-1))) );
%     end
% end
%poisson=poisson;

% for i=1:(x_size-i_x-1)/(2*2)+1
%     for j=1:(y_size-i_y-1)/(2*2)+1
%         poisson(i,j)=poisson(i*2-1,j*2-1);
%     end
% end

%clear poisson
% size(poisson)
W=poisson;
W(:,:)=1;
ZI=polyfitweighted2(Y1,X1,poisson,s,W);
%ZI=polyfitweighted2Mod(Y1,X1,poisson/2,s,smooth_vec);
smooth_comp=polyval2(ZI,Y1,X1);

[m,n]=size(poisson);
%smooth_comp=interp2(1:m,1:n,smooth_comp_first,XI3,YI3,'*spline');


%[alpha_x_DM,alpha_y_DM] = cluster_lense_model_loops_DM;

%display('started fourier DM')

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
% [m,n]=size(poisson);
% alpha_x_DM=interp2(1:m,1:n,alpha_y_DMfs,XI,YI,'*spline');
% alpha_y_DM=interp2(1:m,1:n,alpha_x_DMfs,XI,YI,'*spline');
alpha_x_DM=alpha_x_DMfs*resto^2/pi;
alpha_y_DM=alpha_y_DMfs*resto^2/pi;
%save 'alpha_xtry.txt' alpha_x_DM -ascii
clear  alpha_x_DMfs alpha_y_DMfs
%display('finished DM fourier')

elseif strcmp(SplineOrGaussian,'Gaussian')==1
smooth_comp=poisson;

%display('started fourier DM')
%x_factor(1:2*length(smooth_comp(:,1)),1:2*length(smooth_comp(1,:)))=0;
%y_factor(1:2*length(smooth_comp(:,1)),1:2*length(smooth_comp(1,:)))=0;
%g_factor(1:2*length(smooth_comp(:,1)),1:2*length(smooth_comp(1,:)))=0;
%[Y,X]=meshgrid(-length(smooth_comp(:,1)):length(smooth_comp(:,1))-1,-length(smooth_comp(1,:)):length(smooth_comp(1,:))-1);
%x_factor=X./(X.^2+Y.^2);
% y_factor=Y./(X.^2+Y.^2);
% g_factor=exp(-(X.^2+Y.^2)./(2.0*s*s))./(2.0*pi*s*s);
% x_factor_shift=fftshift(x_factor);
% y_factor_shift=fftshift(y_factor);
% g_factor_shift=fftshift(g_factor);
% x_factor_shift(1,1)=0;
% y_factor_shift(1,1)=0;
% xfourier=fft2(x_factor_shift);
% yfourier=fft2(y_factor_shift);
% gfourier=fft2(g_factor_shift);
%median(median(smooth_comp(1:5,:)))+median(median(smooth_comp(:,1:5)))+median(median(smooth_comp(leng2-5:leng2,:)))+median(median(smooth_comp(:,leng2-5:leng2)))/4
new_smooth(1:2*length(smooth_comp(:,1)), 1:2*length(smooth_comp(1,:)))=0;%median(median(smooth_comp(1:5,:)))+median(median(smooth_comp(:,1:5)))+median(median(smooth_comp(leng2-5:leng2,:)))+median(median(smooth_comp(:,leng2-5:leng2)))/4;
new_smooth(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)))=smooth_comp;

fourier_smooth_new=fft2(new_smooth);
s=params_data(2);
  g_factor=exp(-(X.^2+Y.^2)./(2.0*s*s))./(2.0*pi*s*s);
  g_factor_shift=fftshift(g_factor);
  gfourier=fft2(g_factor_shift);
alpha_x_DMf=ifft2(xfourier.*gfourier.*fourier_smooth_new);
alpha_y_DMf=ifft2(yfourier.*gfourier.*fourier_smooth_new);
%clear new_smooth xfourier yfourier
alpha_x_DMfs=alpha_x_DMf(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)));
alpha_y_DMfs=alpha_y_DMf(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)));
clear alpha_x_DMf alpha_y_DMf
% [XI, YI]=meshgrid(1:0.2:(x_size-i_x-1)/5+1);
% [m,n]=size(poisson);
% alpha_x_DM=interp2(1:m,1:n,alpha_y_DMfs,XI,YI,'*spline');
% alpha_y_DM=interp2(1:m,1:n,alpha_x_DMfs,XI,YI,'*spline');
alpha_x_DM=alpha_x_DMfs*resto^2/pi;
alpha_y_DM=alpha_y_DMfs*resto^2/pi;
%save 'alpha_xtry.txt' alpha_x_DM -ascii
clear  alpha_x_DMfs alpha_y_DMfs
%display('finished DM fourier')

end
%[magnification_ALL]=calculate_magnification_ALL_general(k_new,k_gal,gamma,phi);
mean(mean(abs(alpha_x)));
mean(mean(abs(alpha_x_DM)));
rat=200/mean(mean(abs(alpha_x_DM)));
rat_gal=200/mean(mean(abs(alpha_x)));% alpha_x=alpha_x/max(max(abs(alpha_x)));
% alpha_x_DM=alpha_x_DM/max(max(abs(alpha_x_DM)));
%rat=1/(mean(mean(abs(alpha_x)+abs(alpha_x_DM)))/2)*k_new
factor=k_new;
%rat=200/mean(mean(abs(alpha_x_DM)))*k_new;
gammacos2phi=gamma*cos(2*phi);
gammasin2phi=gamma*sin(2*phi);
alpha_x_add=gammacos2phi*XI2+gammasin2phi*YI2;
alpha_y_add=gammasin2phi*XI2-gammacos2phi*YI2;
%save tillnow1
clear alpha_x_ALL1 alpha_y_ALL1
alpha_x_ALL1=(k_gal*rat_gal*alpha_x(1:leng2,1:leng2)+(1-k_gal)*rat*alpha_x_DM(1:leng2,1:leng2))*k_new+alpha_x_add(1:leng2,1:leng2);%+k_gal*k_new*k_plus*alpha_x_plus(1:leng,1:leng);
alpha_y_ALL1=(k_gal*rat_gal*alpha_y(1:leng2,1:leng2)+(1-k_gal)*rat*alpha_y_DM(1:leng2,1:leng2))*k_new+alpha_y_add(1:leng2,1:leng2);%+k_gal*k_new*k_plus*alpha_y_plus(1:leng,1:leng);
 %alpha_x_ALL1=interp2(alpha_x_ALL,1:0.25:leng2+0.75,(1:0.25:leng2+0.75)');
 %alpha_y_ALL1=interp2(alpha_y_ALL,1:0.25:leng2+0.75,(1:0.25:leng2+0.75)');
% clear alpha_x_ALL1 alpha_y_ALL1
%this factor should be itterated. Always with plus sign: 
%display('finished calculating total def field')
% %calculate magnification and mass:

if strcmp(GPU_OR_CPU,'GPU')==1
alpha_x_ALL=gather(alpha_x_ALL1);
alpha_y_ALL=gather(alpha_y_ALL1);
elseif strcmp(GPU_OR_CPU,'CPU')==1
alpha_x_ALL=(alpha_x_ALL1);
alpha_y_ALL=(alpha_y_ALL1);
end
  [da_x_dy,da_x_dx]=gradient(alpha_x_ALL/resto);
  [da_y_dy,da_y_dx]=gradient(alpha_y_ALL/resto);
  poisson_ALL=da_x_dx+da_y_dy;
  magnification_ALL=abs(1./(1-poisson_ALL+da_x_dx.*da_y_dy-da_x_dy.*da_y_dx));
% % %magnification_ALL_with_sign=(1./(1-poisson_ALL+da_x_dx.*da_y_dy-da_x_dy.*da_y_dx));
  image(fliplr((magnification_ALL))')
%   api=find(magnification_ALL>5);
%   length(api)*0.05^2/3600*resto^2
  waitforbuttonpress
% v=[0.2:0.2:3]
% contour(poisson_ALL,v)
     save('bestModel.mat','alpha_x_ALL','alpha_y_ALL','i_x','i_y','x_start','y_start','leng','x_size','y_size','magnification_ALL');
%%save bestmodel_avgsrc10000_1
%display('finished calculating tmagnification and kappa')
%display('delensing-relensing and calculating chi^2')
%now calculatechi2
%tic


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATING CHI2 AND RMS - NO NEED TO MEDDLE HERE BUT MARK IF 
%YOU WANT TO CHECK PARITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma_location=0.5/pix_scale;
min_disto(1:length(f(:,1)))=10000;
count=0;
flag_lensingok=1;
flag_lensingok2=1;
count=count+1;
no_systems=f(1,5);
if strcmp(GPU_OR_CPU,'GPU')==1
deflection = gpuArray(complex(alpha_x_ALL,alpha_y_ALL));
elseif strcmp(GPU_OR_CPU,'CPU')==1
deflection = (complex(alpha_x_ALL,alpha_y_ALL));    
end

line_no=1;
%sum_disto(1:no_systems)=0;
for i=1:no_systems
    if find(Ind_Free_Reds(:)==i)>0
      atemptu=find(Ind_Free_Reds(:)==i);
       %params_data(atempto+6)7+No_of_Bright_Free+1:Nparams
       if CosmoFree==0
      te=params_data(atemptu+12+No_of_Bright_Free);%%need to add weights
       elseif CosmoFree==1
       te=params_data(atemptu+17+No_of_Bright_Free);%%need to add weights 
       %here the cosmo te is in redshift not dls/ds!!
       end
    else
       te=0;
    end 
% now working on first system
no_images_in_sys=f(line_no,4);
%s_check(1:x_size-i_x,1:y_size-i_y)=0;
countr=0;
for no_images=1:no_images_in_sys
    if flag_lensingok==1  
        countr=countr+1;
%     omegaM
%     omegaL
%     w_0
%     wa
if CosmoFree==1
   [k_systemte]=calculate_k_for_system3(cluster_redshift, main_source_z, (f(line_no,7)+te), H_zero ,omegaM,omegaL,w_0,wa);
else
  k_systemte=f(line_no,3)+te; 
end
        a=((f(line_no,1)-x_start+1));
        b=((f(line_no,2)-y_start+1));
        a1=round((f(line_no,1)-x_start+1)/resto);
        b1=round((f(line_no,2)-y_start+1)/resto);
           magnification_ALL(a1-5:a1+5,b1-5:b1+5)=100;
if find(Ind_Parity_Force(:)==i)>0   
    %if i==1% || i==5 || i==6 %  || i==7  || i==8 || i==15% || i==15 
[da_x_dy,da_x_dx]=gradient(alpha_x_ALL*k_systemte/resto);
[da_y_dy,da_y_dx]=gradient(alpha_y_ALL*k_systemte/resto);
poisson_ALLtemp=da_x_dx+da_y_dy;
%magnification_ALL_temp=abs(1./(1-poisson_ALL+da_x_dx.*da_y_dy-da_x_dy.*da_y_dx));
magnification_ALL_with_sign_temp=(1./(1-poisson_ALLtemp+da_x_dx.*da_y_dy-da_x_dy.*da_y_dx));
if sign(magnification_ALL_with_sign_temp(a1,b1))==f(line_no,6) %&& sign((1-poisson_ALL(a,b)/2))==f(line_no,7)
%then nothing
%display('parity signs are ok');
else
   flag_lensingok2=0; 
end
  end

        for m=a
           for n=b
                im_x(countr)=m-round(alpha_x_ALL(a1,b1)*k_systemte);
                im_y(countr)=n-round(alpha_y_ALL(a1,b1)*k_systemte);
                if im_x(countr)>0 && im_x(countr)<(x_size-i_x-2) && im_y(countr)>0 && im_y(countr)<(y_size-i_y-2)
               %if im_x(countr)>0 && im_x(countr)<(leng2-2) && im_y(countr)>0 && im_y(countr)<(leng2-2)
                else
                    display('could not lens source - out of FOV; with:')
                    %q,s,k_new,k_gal,gamma,phi
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
    clear source_pos
    source_pos = imxavg + 1i * imyavg;
clear distance distance_gpu
    % should compute the actual distance from the source if abs=sqrt(real^2+img^2) as it should be:
    distance_gpu=abs(grid_image-k_systemte*deflection-source_pos);
    if strcmp(GPU_OR_CPU,'GPU')==1
    distance=gather(distance_gpu);
    elseif strcmp(GPU_OR_CPU,'CPU')==1
    distance=(distance_gpu);
    end
    
    
    %image(distance)
    %here should measure the distance from real image
    %for ja=1:no_images_in_sys
    clear asu asu0 m n im_x im_y
    asu0=find(distance<=resto*2);
    iasu=1;
    for p=1:length(asu0)
        [m,n]=ind2sub(size(distance),asu0(p));
        distance0=distance(m,n);
        addpoint=1;
        for im=-1:1
            for in=-1:1
                if (im==0 && in==0)
                    continue
                end
                im1=m+im;
                in1=n+in;
                if im1<1 || im1>leng2 || in1<1 || in1>leng2
                    continue
                end
                if distance(im1,in1)<distance0
                    addpoint=0;
                    break
                end
            end
            if addpoint==0
                break
            end
        end
        if addpoint==1
            asu(iasu)=asu0(p);
            iasu=iasu+1;
        end
    end
   if iasu>1
        for p=1:length(asu)
        [m,n]=ind2sub(size(distance),asu(p));
        
              %%%%% mark to show images
         magnification_ALL(m-5:m+5,n-5:n+5)=30;
         image(fliplr((magnification_ALL))')
         %%%%%%%%%%%%%%%%%%%%%%%%%
         
         for ja=1:no_images_in_sys
            min_distotemp=sqrt( (f(line_no+ja-no_images,1)-x_start+1-m*resto)^2+(f(line_no+ja-no_images,2)-y_start+1-n*resto)^2);
            if min_disto(line_no+ja-no_images)>min_distotemp
                min_disto(line_no+ja-no_images)=min_distotemp;
                %closest_image_m(line_no+ja-no_images)=m*resto+x_start;
                %closest_image_n(line_no+ja-no_images)=n*resto+y_start;
            end
        end
        end  
        else
        flag_lensingok=0;
    end
end

line_no=line_no+1;

end %going through systems

if flag_lensingok==1 && flag_lensingok2==1
min_disto=min_disto*pix_scale;
rms=sqrt(sum((min_disto).^2)/length(f(:,1)));
final_chi2=sum((min_disto).^2/0.5^2);
display(['chi2 is: ' num2str(final_chi2) ', rms ("): ' num2str(rms)])
display('parity OK')
display('-----------------------')
elseif flag_lensingok==1 && flag_lensingok2~=1
min_disto=min_disto*pix_scale;
rms=sqrt(sum((min_disto).^2)/length(f(:,1)))
final_chi2=sum((min_disto).^2/0.5^2)+50;
final_chi2
 else
    final_chi2=10^9;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINISHED CALCULATING CHI2 AND RMS for SL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if relat_weight_WL>0

%%%%%%%%% now check chi2 of fit to WL shape measurements

if CosmoFree==0
[kmeanred]=calculate_k_for_system3(cluster_redshift, main_source_z,mean_sample_redshift ,70,0.3,0.7,-1,0);   %%ATTENZIONE: what is this?
elseif CosmoFree==1
 [kmeanred]=calculate_k_for_system3(cluster_redshift, main_source_z, mean_sample_redshift, H_zero ,omegaM,omegaL,w_0,wa);
end
[da_x_dy,da_x_dx]=gradient(alpha_x_ALL*kmeanred/resto);
[da_y_dy,da_y_dx]=gradient(alpha_y_ALL*kmeanred/resto);
poisson_ALL=da_x_dx+da_y_dy;
magnification_ALL=(1./(1-poisson_ALL+da_x_dx.*da_y_dy-da_x_dy.*da_y_dx));


chi2_WL=0;

g1(1:leng2,1:leng2)=0;
g2(1:leng2,1:leng2)=0;

gamma1=0.5*(da_x_dx-da_y_dy);
gamma2=da_y_dx;
gamma_abs=sqrt(gamma1.^2+gamma2.^2);
A=(gamma_abs./(1-poisson_ALL/2)<1);
B=(gamma_abs./(1-poisson_ALL/2)>=1);
g1(A)=gamma1(A)./(1-poisson_ALL(A)/2); % this is g_1
g1(B)=gamma1(B).*(1-poisson_ALL(B)/2)./gamma_abs(B).^2;

g2(A)=gamma2(A)./(1-poisson_ALL(A)/2); % this is g_2
g2(B)=gamma2(B).*(1-poisson_ALL(B)/2)./gamma_abs(B).^2;

f_length=1; %in pixels, fiudducial length for g normalization; %ignore

g1cat(1:leng2,1:leng2)=0;
g2cat(1:leng2,1:leng2)=0; 

%for plotting:
gxcat(1:leng2,1:leng2)=0;
gycat(1:leng2,1:leng2)=0;
x2cat(1:leng2,1:leng2)=0;
y2cat(1:leng2,1:leng2)=0;

chi2_WL=0;
flag_entered=0;
countthrown=0;
accepted=0;

for t=1:length(fwl(:,1))

    xpos=round( (fwl(t,1)/pix_scale+bcg_x)/resto); %in pixels in shortended LR frame (LR by restoxresto as usual)
    ypos=round( (fwl(t,2)/pix_scale+bcg_y)/resto); %in pixels in shortended LR frame (LR by restoxresto as usual)
    

    if xpos>0 && xpos<leng2 && ypos>0 && ypos<leng2
        
        g1cat(xpos,ypos)=fwl(t,5);  %observed ellipticity plus component; this is e1
        g2cat(xpos,ypos)=fwl(t,6); %observed ellipticity cross component ; e2
        x2cat(xpos,ypos)=xpos;
        y2cat(xpos,ypos)=ypos;
        SNRgal(t)=fwl(t,4);
        sigma_plus(t)=0.32;%sqrt(0.5*(0.3^2+(3/SNRgal(t))^2));%%%%0.3/sqrt(mean_num_obj_pixel);
        sigma_cross(t)=sigma_plus(t);%%%0.3/sqrt(mean_num_obj_pixel);

      %calculate chi2 WL:
         if  magnification_ALL(xpos,ypos)>0 && (poisson_ALL(xpos,ypos)/2)<1 
         chi2_WL= chi2_WL+( (g1cat(xpos,ypos)-g1(xpos,ypos)).^2 ./sigma_plus(t).^2 + (g2cat(xpos,ypos)-g2(xpos,ypos)).^2 ./sigma_cross(t).^2);%*(1-(poisson_ALL(xpos,ypos)/2)^2) ;
         accepted=accepted+1;
         flag_entered=1;
         elseif magnification_ALL(xpos,ypos)<0 || (poisson_ALL(xpos,ypos)/2 )>=1  
              countthrown=countthrown+1;
        
          end  
      
   
    end
end


num_images=length(f(:,1)); %num of multiple images
num_systems=f(1,5); %num of systems

num_pixels_shear=round(leng2/jumpi)^2;

ndf_SL=2*(num_images-num_systems)-num_param; %degrees of freedom SL
ndf_WL=2*accepted-num_param; %degrees of freedom WL

if DOFnormalize==1
if flag_entered==1 && final_chi2<10^9

    chi2red_WL=chi2_WL/ndf_WL;
    chi2red_SL=final_chi2/ndf_SL;
       
    final_chi2=fiducial_chi2_mult*(relat_weight_SL*(final_chi2/ndf_SL)+relat_weight_WL*chi2_WL/ndf_WL);%/(10^7)  %%ATTENZIONE parametro 100*!
    final_chi2

elseif flag_entered==0 || final_chi2>10^8.8
    
    chi2_WL=10^9;
      final_chi2=fiducial_chi2_mult*(relat_weight_SL*(final_chi2/ndf_SL)+relat_weight_WL*chi2_WL/ndf_WL);%/(10^7)  %%ATTENZIONE parametro 100*!
  final_chi2
    
end
else
 if flag_entered==1 && final_chi2<10^9

    chi2red_WL=chi2_WL/ndf_WL;
    chi2red_SL=final_chi2/ndf_SL;
       
    final_chi2=fiducial_chi2_mult*(relat_weight_SL*(final_chi2)+relat_weight_WL*chi2_WL);%/(10^7)  %%ATTENZIONE parametro 100*!
    final_chi2

elseif flag_entered==0 || final_chi2>10^8.8
    
    chi2_WL=10^9;
      final_chi2=fiducial_chi2_mult*(relat_weight_SL*(final_chi2)+relat_weight_WL*chi2_WL);%/(10^7)  %%ATTENZIONE parametro 100*!
  final_chi2
    
end   
end
else
num_images=length(f(:,1)); %num of multiple images
num_systems=f(1,5); %num of systems

num_pixels_shear=round(leng2/jumpi)^2;

ndf_SL=2*(num_images-num_systems)-num_param; %degrees of freedom SL
ndf_WL=2*accepted-num_param; %degrees of freedom WL
 
    
end

%%%% now if you wish to plot the shera and average unmark the following:
% %%%ATTENZIONE: begin added by AGNESE
%accepted=0;
%chi2_WL=0;
% %%%%%%%%% now check chi2 of fit to WL shape measurements
% % remind that:
% %define temporarily: 
% %kmeanred=1;
% %mean_sample_redshift=1.12;  %%ATTENZIONE parametro!
% if CosmoFree==0
% [kmeanred]=calculate_k_for_system3(cluster_redshift, main_source_z,mean_sample_redshift ,70,0.3,0.7,-1,0);   %%ATTENZIONE: what is this?
% elseif CosmoFree==1
%  [kmeanred]=calculate_k_for_system3(cluster_redshift, main_source_z, mean_sample_redshift, H_zero ,omegaM,omegaL,w_0,wa);
% end   
% [da_x_dy,da_x_dx]=gradient(alpha_x_ALL*kmeanred/resto);
% [da_y_dy,da_y_dx]=gradient(alpha_y_ALL*kmeanred/resto);
% poisson_ALL=da_x_dx+da_y_dy;
% magnification_ALL=(1./(1-poisson_ALL+da_x_dx.*da_y_dy-da_x_dy.*da_y_dx));
% 
% % assume shape catalog as follows: 
% % col 1 = index ; 
% % col 2 = x (in pixels for now);
% % col 3 = y (in pixels for now);
% % col 4= shear_plus; 
% % col 5 = shear_cross
% %f=load(shear_catalog);
% 
% chi2_WL=0;
% 
% %fwl=load('JulianCatMS2137shapes.txt');  %%ATTENZIONE: nome catalogo da usare
% 
% %for now assuming constant redshift of background: z=1.5 but need to test
% %each time or take for each galaxy its redshift!
% %need to check resto factors are correct
% elk_plus(1:leng2,1:leng2)=0;
% elk_cross(1:leng2,1:leng2)=0;
% 
% gamma_plus=0.5*(da_x_dx-da_y_dy);
% gamma_cross=da_y_dx;
% gamma_abs=sqrt(gamma_plus.^2+gamma_cross.^2);
% A=(gamma_abs./(1-poisson_ALL/2)<1);
% B=(gamma_abs./(1-poisson_ALL/2)>=1);
% elk_plus(A)=gamma_plus(A)./(1-poisson_ALL(A)/2);
% elk_plus(B)=gamma_plus(B).*(1-poisson_ALL(B)/2)./gamma_abs(B).^2;
% 
% elk_cross(A)=gamma_cross(A)./(1-poisson_ALL(A)/2);
% elk_cross(B)=gamma_cross(B).*(1-poisson_ALL(B)/2)./gamma_abs(B).^2;
% 
% 
% %jumpi=60*10/resto; %salto  %%%ATTENZIONE parametro!!!
% 
% c_start=round(jumpi/2)+1;
% c_end=leng2-round(jumpi/2)+1;
% 
% for i=c_start:jumpi:c_end%th(f(:,1)) 
%     for j=c_start:jumpi:c_end
%         x2(i,j)=XI_ram2(i,j);%round((bcg_x+f(i,1)/pix_scale)/resto);  %%%ATTENZIONE: è giusto??
%         y2(i,j)=YI_ram2(i,j);%round((bcg_y+f(i,2)/pix_scale)/resto);
% 
%         % consider averaging also for model
%     
%         G_plus(i,j)=elk_plus(i,j);  %reduced shear plus component
%         G_cross(i,j)=elk_cross(i,j);  %reduced shear cross component
%         
%         %and don't forget to account for different background redshift!!!
%         G_abs(i,j)=sqrt(G_plus(i,j).^2+G_cross(i,j).^2);
%         phi_abs(i,j)=0.5*atan2(G_cross(i,j),G_plus(i,j));
%         
%         tempx(i,j)=G_abs(i,j).*cos(phi_abs(i,j)); %reduced shear x-component
%         tempy(i,j)=G_abs(i,j).*sin(phi_abs(i,j)); %reduced shear y-component
%     end
% end
% 
% %now compare to shear catalog:
% %sigma_plus(i)=0.3*f(i,3);
% %sigma_cross(i)=0.3*f(i,4);
% %chi2_WL= chi2_WL+(f(i,3)-gamma_plus(i))^2/sigma_plus(i)^2 + (f(i,4)-gamma_cross(i))^2/sigma_cross(i)^2 ;
% 
% f_length=1; %in pixels, fiudducial length for g normalization
% 
% %and don't forget to account for different background redshift!!!
% % display('now saving')
% %waitforbuttonpress
% %close
% %figure
% % %streamslice(51:100:leng-49,51:100:leng-49,(tempx(51:100:leng-49,51:100:leng-49)),(tempy(51:100:leng-49,51:100:leng-49)),6);
% 
% 
% 
% quiver( x2(c_start:jumpi:c_end, c_start:jumpi:c_end), y2(c_start:jumpi:c_end, c_start:jumpi:c_end), tempx(c_start:jumpi:c_end, c_start:jumpi:c_end), tempy(c_start:jumpi:c_end, c_start:jumpi:c_end),1,'blue');
% %waitforbuttonpress
% hold on
% 
% G_pluscat(1:leng2,1:leng2)=0;
% G_crosscat(1:leng2,1:leng2)=0; 
% G_abscat(1:leng2,1:leng2)=0;
% phi_abscat(1:leng2,1:leng2)=0;
% 
% tempxcat(1:leng2,1:leng2)=0;
% tempycat(1:leng2,1:leng2)=0;
% x2cat(1:leng2,1:leng2)=0;
% y2cat(1:leng2,1:leng2)=0;
% 
% for t=1:length(fwl(:,1))
% 
%     xpos=round( (fwl(t,1)/pix_scale+bcg_x)/resto);
%     ypos=round( (fwl(t,2)/pix_scale+bcg_y)/resto);
%     
% 
%     if xpos>0 && xpos<leng2 && ypos>0 && ypos<leng2
%         
%         G_pluscat(xpos,ypos)=fwl(t,5);  %observed ellipticity plus component
%         G_crosscat(xpos,ypos)=fwl(t,6); %observed ellipticity cross component
%         x2cat(xpos,ypos)=xpos;
%         y2cat(xpos,ypos)=ypos;
%         
%         %and don't forget to account for different background redshift!!!
%         G_abscat(xpos,ypos)=sqrt(G_pluscat(xpos,ypos).^2+G_crosscat(xpos,ypos).^2); 
%         phi_abscat(xpos,ypos)=0.5*atan2(G_crosscat(xpos,ypos),G_pluscat(xpos,ypos));
%        
%         tempxcat(xpos,ypos)=G_abscat(xpos,ypos).*cos(phi_abscat(xpos,ypos));  %observed ell x-component
%         tempycat(xpos,ypos)=G_abscat(xpos,ypos).*sin(phi_abscat(xpos,ypos));  %observed ell y-component
%         
%     end
% end
% 
% 
% 
% %%% now average on big pixels of size jumpi
% 
% xposthrow(1:leng2,1:leng2)=0;
% yposthrow(1:leng2,1:leng2)=0;
% tempxthrow(1:leng2,1:leng2)=0;
% tempythrow(1:leng2,1:leng2)=0;
% 
% x21(1:leng2,1:leng2)=0;
% y21(1:leng2,1:leng2)=0;
%     
% num_pixels_shear=round(leng2/jumpi)^2;
% mean_num_obj_pixel=length(fwl(:,1))/num_pixels_shear;
% 
% sigma_plus=0.3/sqrt(mean_num_obj_pixel);
% sigma_cross=0.3/sqrt(mean_num_obj_pixel);
% 
% 
% chi2_WL=0;
% flag_entered=0;
% countthrown=0;
% accepted=0;
%  
% 
% conta=0;
% kr=0;
% 
% for i=c_start:jumpi:c_end %%th(f(:,1)) 
%  
%     kr=kr+1;
%     kc=0;
%       
%       
%     for j=c_start:jumpi:c_end 
% 
%           conta=conta+1;   
%           kc=kc+1;
%     
%           x21(i,j)=XI_ram2(i,j);%round((bcg_x+f(i,1)/pix_scale)/resto);
%           y21(i,j)=YI_ram2(i,j);%round((bcg_y+f(i,2)/pix_scale)/resto);
%           xpos=XI_ram2(i,j);
%           ypos=YI_ram2(i,j);
%    
%           
%           m=G_pluscat(i-round(jumpi/2):i+round(jumpi/2)-1,j-round(jumpi/2):j+round(jumpi/2)-1); %faccio media con comp cross e plus
%           g_mean_plus(i,j)=sum(sum(m))./sum(sum(m~=0));
%           nobj_plus(kr,kc)=sum(sum(m~=0));
%                  
%           m=G_crosscat(i-round(jumpi/2):i+round(jumpi/2)-1,j-round(jumpi/2):j+round(jumpi/2)-1);  %faccio media con comp cross e plus
%           g_mean_cross(i,j)=sum(sum(m))./sum(sum(m~=0));
%           nobj_cross(kr,kc)=sum(sum(m~=0));
%            
%            
%           g_mean_abs(i,j)=sqrt(g_mean_plus(i,j).^2+g_mean_cross(i,j).^2);
%           phi_mean(i,j)=0.5*atan2(g_mean_cross(i,j),g_mean_plus(i,j));
%           
%           tempx1(i,j)=g_mean_abs(i,j)*cos(phi_mean(i,j));  
%           tempy1(i,j)=g_mean_abs(i,j)*sin(phi_mean(i,j));
% 
%          
%           
%           if  magnification_ALL(xpos,ypos)>0 && (poisson_ALL(xpos,ypos)/2)<1 && (nobj_plus(kr,kc)+nobj_cross(kr,kc))~=0 
%               
%               residuo_x(conta)=tempx1(xpos,ypos)-tempx(xpos,ypos);
%               residuo_y(conta)=tempy1(xpos,ypos)-tempy(xpos,ypos);
%               residuo_chi(conta)=sqrt( (tempx1(xpos,ypos)-tempx(xpos,ypos)).^2/sigma_plus^2 + (tempy1(xpos,ypos)-tempy(xpos,ypos)).^2/sigma_cross^2 )*(1-(poisson_ALL(xpos,ypos)/2)^2);
%               flag_entered=1;
%               
%               %%%ATTENZIONE: quali componenti voglio usare nel chi2_WL????
%               %chi2_WL= chi2_WL+( (G_plus1(xpos,ypos)-G_plus(xpos,ypos)).^2 ./sigma_plus.^2 + (G_cross1(xpos,ypos)-G_cross(xpos,ypos)).^2 ./sigma_cross.^2)*(1-(poisson_ALL(xpos,ypos)/2)^2) ;
%               chi2_WL= chi2_WL+( (tempx1(xpos,ypos)-tempx(xpos,ypos)).^2 ./sigma_plus.^2 + (tempy1(xpos,ypos)-tempy(xpos,ypos)).^2 ./sigma_cross.^2)*(1-(poisson_ALL(xpos,ypos)/2)^2) ;
%               accepted=accepted+1;
%           
%           elseif magnification_ALL(xpos,ypos)<0 || (poisson_ALL(xpos,ypos)/2 )>=1 || (nobj_plus(kr,kc)+nobj_cross(kr,kc))==0 
%           
%               countthrown=countthrown+1;
%               xposthrow(xpos,ypos)=xpos;
%               yposthrow(xpos,ypos)=ypos;
%               tempxthrow(xpos,ypos)=tempx(xpos,ypos);
%               tempythrow(xpos,ypos)=tempy(xpos,ypos);
%               
%           end
%         
%      end
% end
% 
% display('chi2 only WL:')
% chi2_WL
%   
% 
% display('number big pixels:')
% conta
% 
% countthrown
% accepted
% 
% quiver( x21(c_start:jumpi:c_end, c_start:jumpi:c_end), y21(c_start:jumpi:c_end, c_start:jumpi:c_end), tempx1(c_start:jumpi:c_end, c_start:jumpi:c_end), tempy1(c_start:jumpi:c_end, c_start:jumpi:c_end),1,'red');
% 
% quiver(xposthrow,yposthrow,tempxthrow,tempythrow,50,'green');
% hold off
% 
% 
% 
% 
% % %%%%agnese test
% % waitforbuttonpress
% % 
% % 
% % [nelements,xcenters]=hist(residuo_x(:),20)
% % hist(residuo_x(:),20)
% % axis([-1 1 0 20])
% % std(residuo_x(:))
% % hi = findobj(gca,'Type','patch');
% % set(hi,'FaceColor','m','EdgeColor','b')
% % clear hi
% 
% % waitforbuttonpress
% % [nelements, xcenters]=hist(residuo_y(:),20)
% % hist(residuo_y(:),20)
% % axis([-1 1 0 10])
% % std(residuo_y(:))
% % hii = findobj(gca,'Type','patch');
% % set(hii,'FaceColor','y','EdgeColor','b')
% % 
% % 
% % lbin=0.5
% % xo=0.0
% % xcent(1)=xo
% % for i=2:1:30
% %    xcent(i)=xcent(i-1)+lbin; 
% % end
% % 
% % waitforbuttonpress
% % [nelements, xxx]=hist(residuo_chi(:),xcent)
% % hist(residuo_chi(:),xcent)
% % axis([-0.1 15 0 20])
% % std(residuo_chi(:))
% % hii = findobj(gca,'Type','patch');
% % set(hii,'FaceColor','r','EdgeColor','b')
% 
% 
% 
% % waitforbuttonpress
% % 
% % scatter(G_plus(:),G_plus1(:), 'black'), hold on
% % plot(min(G_plus(:)):0.01:max(G_plus(:)),min(G_plus(:)):0.01:max(G_plus(:)),'red')
% % scatter(G_cross(:),G_cross1(:),'magenta'), hold off
% % waitforbuttonpress
% % hist(G_plus1(:),30)
% % h = findobj(gca,'Type','patch');
% % set(h,'FaceColor','r','EdgeColor','r')
% % clear h
% % waitforbuttonpress
% % hist(G_plus(:),30)
% 
% 
% num_images=length(f(:,1)); %num of multiple images
% num_systems=f(1,5); %num of systems
% %num_param=6; %num of parameters  %%ATTENZIONE parametro!
% num_pixels_shear=round(leng2/jumpi)^2;
% 
% ndf_SL=2*(num_images-num_systems)-num_param; %degrees of freedom SL
% ndf_WL=2*num_pixels_shear-num_param; %degrees of freedom WL
% %relat_weight_SW=0.02; % weight SL respect to WL in chi2 minimization   %%ATTENZIONE parametro!
% 
% 
% 
% if flag_entered==1 && final_chi2<10^9
% 
%     chi2red_WL=chi2_WL/ndf_WL;
%     chi2red_SL=final_chi2/ndf_SL;
%     
%         
%     final_chi2=fiducial_chi2_mult*(relat_weight_SL*(final_chi2/ndf_SL)+relat_weight_WL*chi2_WL/ndf_WL);%/(10^7)  %%ATTENZIONE parametro 100*!
%     final_chi2
% 
% elseif flag_entered==0 || final_chi2>10^8.8
%     
%     chi2_WL=10^9;
%       final_chi2=fiducial_chi2_mult*(relat_weight_SL*(final_chi2/ndf_SL)+relat_weight_WL*chi2_WL/ndf_WL);%/(10^7)  %%ATTENZIONE parametro 100*!
%   final_chi2
%     
% end
% 
% 
% 
% 
% 
% %%%ATTENZIONE: end added by AGNESE
