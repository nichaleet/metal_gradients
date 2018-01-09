function make_catalog_with_abmag;
redshift = 0.465;
f=load('temp_cat_all.txt');
lines1=length(f(:,1));
% #   1 NUMBER          Running object number
% #   2 X_IMAGE         Object position along x                         [pixel]
% #   3 Y_IMAGE         Object position along y                         [pixel]
% #   4 MAG_ISO         Isophotal magnitude of 475                      [mag]
% #   5 FLUX_ISO        Isophotal flux      of 475                      [count]
% #   6 MAG_ISO         Isophotal magnitude of 814                      [mag]
% #   7 FLUX_ISO        Isophotal flux      of 814                      [count]
% #   8 R-B             mag of 814-mag of 475                           [mag]
% #   9 FWHM_IMAGE      FWHM assuming a gaussian core                   [pixel]
% #  10 ALPHA_J2000     Right ascension of barycenter (J2000)           [deg]
% #  11 DELTA_J2000     Declination of barycenter (J2000)               [deg]
% #  12 A_IMAGE         Profile RMS along major axis                    [pixel]
% #  13 B_IMAGE         Profile RMS along minor axis                    [pixel]
% #  14 THETA_IMAGE     Position angle (CCW/x)                          [deg]
% #  15 ELLIPTICITY     1 - B_IMAGE/A_IMAGE
[dl,dm] = lum_dist(redshift);
 new_file(:,1:6)=f(:,1:6);
 new_file(:,7)=f(:,6)-dm;
 new_file(:,8:15)=f(:,8:15);

save temp_cat_all_withAbmag_z0465.txt new_file -ascii;
