pro magnification
;calculate magnification for cswa19
;the magnification map is from matlab
;>> load('bestModel.mat');
;>> fitswrite(fliplr(imrotate(magnification_ALL,270)),'magnification1.fits');
magmap = readfits('Magnification_CSWA19_2p0325_20140508_LTM_Gc.fits')
detectionmap = readfits('CSWA19handmosaic_Ha_detection_interp.fits')
detectionmap = detectionmap[0:1199,0:1199]

good = where(detectionmap gt 5.)
print,'by area', mean(magmap(good))

;flux weighted
imha = readfits('CSWA19handmosaic_Ha_interp.fits')
imha = imha[0:1199,0:1199]
good = where(detectionmap gt 5. and imha ne 0.)
mag = total(imha(good))/total(imha(good)/magmap(good))
print, 'flux weighted:',mag
stop



end
