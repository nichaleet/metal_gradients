pro magnification
mag =readfits('magnification1.fits')
detect = readfits('CSWA159_Ha_detection_interp.fits')
good = where(detect gt 10.)

print,'By area:', mean(mag(good))

;flux weighted
imha = readfits('CSWA159_Ha_interp.fits')
good = where(detect gt 10. and imha ne 0.)
mag = total(imha(good))/total(imha(good)/mag(good))
print, 'Flux weighted:',mag
stop

end
