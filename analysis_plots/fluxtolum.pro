function fluxtolum,flux,z
;This is for flux in unit of erg/s/cm^2 (not per Hz)
dist = lumdist(z) ;Mpc
dist = double(dist)
lum  = flux*4.*!dpi*(dist*1.d6*3.0857d18)^2
return,lum
end
