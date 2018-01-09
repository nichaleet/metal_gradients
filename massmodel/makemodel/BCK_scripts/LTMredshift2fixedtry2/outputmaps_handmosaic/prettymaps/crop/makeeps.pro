pro makeeps
;cswa19 is at z=2.03 whose scale is 8.461 kpc/" or 0.039618 " per pixel But the image is twice the angular size so the final scale is 0.167 kpc/pixel

outputdir='/scr2/nichal/workspace/paperplots/cswa19'   
;fitstoeps,file,scale,xtitle,ytitle,colorbartitle,directory_out

;Halpha
fitstoeps,'sourceha_pretty.fits',0.167,'kpc','cswa19','Ha flux(erg/s/cm2)',outputdir,format='(F6.3)'
;Halpha_detection
fitstoeps,'sourceha_detection_pretty.fits',0.167,'kpc','cswa19','Ha detection(sigma)',outputdir,format='(F6.3)'
;Nii
fitstoeps,'sourcenii_pretty.fits',0.167,'kpc','cswa19','[NII] flux(erg/s/cm2)',outputdir,format='(F6.3)'
;OIII
fitstoeps,'sourceoiii_pretty.fits',0.167,'kpc','cswa19','[OIII] flux(erg/s/cm2)',outputdir,format='(F6.3)'
;Hbeta
fitstoeps,'sourcehb_pretty.fits',0.167,'kpc','cswa19','Hb flux(erg/s/cm2)',outputdir,format='(F6.3)'
;kinematic
fitstoeps,'sourcekinematic_pretty.fits',0.167,'kpc','cswa19','velocity(km/s)',outputdir,format='(I0)'
;dispersion
fitstoeps,'sourceveldisp_pretty.fits',0.167,'kpc','cswa19','velocity dispersion(km/s)',outputdir,format='(I0)'
;N2index
fitstoeps,'sourcen2index_pretty.fits',0.167,'kpc','cswa19','[NII]/Ha',outputdir,format='(F5.2)'
;Bayesianmetal
;N2index
fitstoeps,'sourceBayesianmetal_pretty.fits',0.167,'kpc','cswa19','12+logO/H(M08)',outputdir,format='(F5.2)'

end
