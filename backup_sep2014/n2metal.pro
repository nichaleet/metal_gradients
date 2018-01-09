pro N2metal,halpha,halpha_err,NII,NII_err,N2metal,N2metal_err,upper_tags

;make a map with NII with upper limit
NII_current = NII
NII_upper = NII_err+NII
change2upper = where(NII lt 0. and NII_upper gt 0.)
NII_current(change2upper) = NII_upper(change2upper)
;stop
upper_tags = halpha*0.
upper_tags(where(finite(upper_tags))) = 1.
upper_tags(change2upper) = 2.
nan = where(NII_upper lt 0.)
upper_tags(nan) = 3.
stop

N2      = ALOG10(NII_current/halpha)
N2_err  = 0.43429*sqrt((NII_err/NII_current)^2+(halpha_err/halpha)^2)
N2metal = 8.90+0.57*N2   ;metallicity unit = 12+log(O/H)
N2metal_err = 0.57*N2_err


end
 
