pro get_bcg,cat,posx,posy
readcol,cat,v1,V2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15, $
            v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28
for i=0,n_elements(posx)-1 do begin
x=posx(i)
y=posy(i)
dis = (v2-x)^2+(v3-y)^2
number =  where(dis eq min(dis))
el  = (v21(number)-v22(number))/(v21(number)+v22(number))
 endfor

;3 Brightest gals
print, '10 Brightest galaxies'
galregion = where(abs(v2-posx(0)) lt 150. and abs(v3-posy(0)) lt 150.)
v8temp = v8(galregion)

bright = []
for ii=0,9 do begin
   bright = [bright,where(v8temp eq min(v8temp))]
   v8temp(bright)= 99.
endfor
print, 'number',v1(galregion(bright))
print, 'x:',v2(galregion(bright))
print, 'y:',v3(galregion(bright))
print, 'mag', v8(galregion(bright))
print, 'elipticity', (v21(galregion(bright))-v22(galregion(bright)))/(v21(galregion(bright))+v22(galregion(bright)))
print, 'angle', v23(galregion(bright))
;stop
end

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
;;#  24 ELLIPTICITY            1 - B_IMAGE/A_IMAGE                                       
;#  25 CLASS_STAR             S/G classifier output                                     
;#  26 FLAGS                  Extraction flags                                          
;#  27 XWIN_IMAGE             Windowed position estimate along x                         [pixel]
;#  28 YWIN_IMAGE             Windowed position estimate along y                         [pixel]
