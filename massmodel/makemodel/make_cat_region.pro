pro make_cat_region,cat,regname

;#   1 NUMBER                 Running object number                                     
;#   2 X_IMAGE                Object position along x                                    [pixel]
;#   3 Y_IMAGE                Object position along y                                    [pixel]
;#   4 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
;#   5 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
;#   6 MAG_AUTO               Kron-like elliptical aperture magnitude                    [mag]
;#   7 MAGERR_AUTO            RMS error for AUTO magnitude                               [mag]
;#   8 MAG_APER               Fixed aperture magnitude vector                            [mag]
;#   9 MAGERR_APER            RMS error vector for fixed aperture mag.                   [mag]
;#  10 MAG_ISO                Isophotal magnitude                                        [mag]
;#  11 MAGERR_ISO             RMS error for isophotal magnitude                          [mag]
;#  12 FLUX_AUTO              Flux within a Kron-like elliptical aperture                [count]
;#  13 FLUXERR_AUTO           RMS error for AUTO flux                                    [count]
;#  14 FLUX_APER              Flux vector within fixed circular aperture(s)              [count]
;#  15 FLUXERR_APER           RMS error vector for aperture flux(es)                     [count]
;#  16 FLUX_ISO               Isophotal flux                                             [count]
;#  17 FLUXERR_ISO            RMS error for isophotal flux                               [count]
;#  18 FWHM_IMAGE             FWHM assuming a gaussian core                              [pixel]
;#  19 FLUX_RADIUS            Fraction-of-light radii                                    [pixel]
;#  20 KRON_RADIUS            Kron apertures in units of A or B                         
;#  21 A_IMAGE                Profile RMS along major axis                               [pixel]
;#  22 B_IMAGE                Profile RMS along minor axis                               [pixel]
;#  23 THETA_IMAGE            Position angle (CCW/x)                                     [deg]
;#  24 ELLIPTICITY            1 - B_IMAGE/A_IMAGE                                       
;#  25 CLASS_STAR             S/G classifier output                                     
;#  26 FLAGS                  Extraction flags                                          
;#  27 XWIN_IMAGE             Windowed position estimate along x                         [pixel]
;#  28 YWIN_IMAGE             Windowed position estimate along y                         [pixel]

readcol,cat,v1,V2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15, $
            v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28
write_ds9_regionfile_xy,v2,v3,filename=regname


end
