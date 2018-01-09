pro magnification
mag =readfits('magnification1.fits')
detect = readfits('CSWA15zoomESI_Ha_detection_interp.fits')
detect = detect[1:825,1:825]
good = where(detect gt 10.)
print,'By area:', mean(mag(good))

;flux weighted
imha = readfits('CSWA15zoomESI_Ha_interp.fits')
imha = imha[1:825,1:825]
good = where(detect gt 5. and imha ne 0.)
mag = total(imha(good))/total(imha(good)/mag(good))
print, 'Flux weighted:',mag
stop


end
