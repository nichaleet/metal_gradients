pro magnification
mag =readfits('magnification1.fits')
mag = mag[1:730,1:730]
detect = readfits('CSWA11_secondgal_Ha_detection_interp.fits')
good = where(detect gt 5.)
print,'By area:', mean(mag(good))

;flux weighted
imha = readfits('CSWA11_secondgal_Ha_interp.fits')
good = where(detect gt 5. and imha ne 0.)
mag = total(imha(good))/total(imha(good)/mag(good))
print, 'Flux weighted:',mag
stop

end
