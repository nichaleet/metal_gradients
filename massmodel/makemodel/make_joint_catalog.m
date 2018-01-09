function make_joint_catalog;
f=load('cswa15_B.cat');
g=load('cswa15_I.cat');
lines1=length(f(:,1));
lines2=length(g(:,1));
% #   1 NUMBER          Running object number
% #   2 X_IMAGE         Object position along x                         [pixel]
% #   3 Y_IMAGE         Object position along y                         [pixel]
% #   4 ALPHA_J2000     Right ascension of barycenter (J2000)           [deg]
% #   5 DELTA_J2000     Declination of barycenter (J2000)               [deg]
% #   6 MAG_AUTO        Kron-like elliptical aperture magnitude         [mag]
% #   7 MAGERR_AUTO     RMS error for AUTO magnitude                    [mag]
% #   8 MAG_APER        Fixed aperture magnitude vector                 [mag]
% #   9 MAGERR_APER     RMS error vector for fixed aperture mag.        [mag]
% #  10 MAG_ISO         Isophotal magnitude                             [mag]
% #  11 MAGERR_ISO      RMS error for isophotal magnitude               [mag]
% #  12 FLUX_AUTO       Flux within a Kron-like elliptical aperture     [count]
% #  13 FLUXERR_AUTO    RMS error for AUTO flux                         [count]
% #  14 FLUX_APER       Flux vector within fixed circular aperture(s)   [count]
% #  15 FLUXERR_APER    RMS error vector for aperture flux(es)          [count]
% #  16 FLUX_ISO        Isophotal flux                                  [count]
% #  17 FLUXERR_ISO     RMS error for isophotal flux                    [count]
% #  18 FWHM_IMAGE      FWHM assuming a gaussian core                   [pixel]
% #  19 FLUX_RADIUS     Fraction-of-light radii                         [pixel]
% #  20 KRON_RADIUS     Kron apertures in units of A or B
% #  21 A_IMAGE         Profile RMS along major axis                    [pixel]
% #  22 B_IMAGE         Profile RMS along minor axis                    [pixel]
% #  23 THETA_IMAGE     Position angle (CCW/x)                          [deg]
% #  24 ELLIPTICITY     1 - B_IMAGE/A_IMAGE
% #  25 CLASS_STAR      S/G classifier output
% #  26 FLAGS           Extraction flags
% #  27 XWIN_IMAGE      Windowed position estimate along x              [pixel]
% #  28 YWIN_IMAGE      Windowed position estimate along y 
%for i=1:lines1
   %for j=1:lines2
    %if ((f(i,1)==g(j,1)) &&(f(i,2)==g(j,2))) %&& (file_data(obj,6)>=15.5) && (file_data(obj,6)<=23))
 new_file(:,1:3)=f(:,1:3);
 new_file(:,4)=f(:,10);
 new_file(:,5)=f(:,16);
 new_file(:,6)=g(:,10);
 new_file(:,7)=g(:,16);
 new_file(:,8)=new_file(:,4)-new_file(:,6);
 new_file(:,9)=g(:,18);
 new_file(:,10:11)=g(:,4:5);
 new_file(:,12:15)=g(:,21:24);
    %end
   %end
   %end
save temp_cat_all.txt new_file -ascii;
%file_data=load('model_colors.txt');
%scatter(file_data(:,7),file_data(:,8));