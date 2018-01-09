pro magnification
mag =readfits('magnification1.fits')
detect = readfits('CSWA139zoomESI_Ha_detection_interp.fits')
detect = detect[0:462,0:462]
good = where(detect gt 10.)

print,'By area:', mean(mag(good))

;flux weighted
imha = readfits('CSWA139zoomESI_Ha_interp.fits')
imha = imha[0:462,0:462]
good = where(detect gt 10. and imha ne 0.)
mag = total(imha(good))/total(imha(good)/mag(good))
print, 'Flux weighted:',mag
stop

end
