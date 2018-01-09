pro magnification
;calculate magnification for cswa19
magmap = readfits('Magnification_CSWA19_2p0325_20140508_LTM_Gc.fits')
detectionmap = readfits('CSWA19handmosaic_Ha_detection_interp.fits')
detectionmap = detectionmap[0:1199,0:1199]

good = where(detectionmap gt 10.)
mag = magmap(good)
print, summary(mag)
print, mean(mag)
stop

end
