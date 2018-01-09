pro O3N2metal,OIII,hbeta,NII,halpha,OIII_err,hbeta_err,NII_err,halpha_err,O3N2metal,O3N2metal_err

;find detected NII pixels
NII_current = NII
NII_upper = NII_err+NII
change2upper = where(NII lt 0. and NII_upper gt 0.)
if change2upper(0) ne -1 then NII_current(change2upper) = NII_upper(change2upper)
nan = where(NII_upper lt 0.)
if nan(0) ne -1 then NII_current(nan) = 1./0.

;find detected Hbeta pixels
Hbeta_current = Hbeta
Hbeta_upper = Hbeta_err+Hbeta
change2upper = where(Hbeta lt 0. and Hbeta_upper gt 0.)
if change2upper(0) ne -1 then Hbeta_current(change2upper) = Hbeta_upper(change2upper)
nan = where(Hbeta_upper lt 0.)
if nan(0) ne -1 then Hbeta_current(nan) = 1./0.

if where(OIII lt 0. and NII lt 0.) ne -1 then stop ; This is to prevent when there are 2 negatives (non-detections) combine to give a false metallicity

;Now, real O3N2 calculation

O3N2 = alog10((OIII/hbeta_current)/(NII_current/halpha))
O3N2metal = 8.73-0.32*O3N2

O3N2_err = 0.43429*sqrt((OIII_err/OIII)^2+(hbeta_err/hbeta_current)^2+(NII_err/NII_current)^2+(halpha_err/halpha)^2)
O3N2metal_err = 0.32*O3N2_err

end
