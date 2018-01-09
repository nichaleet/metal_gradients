pro makeeps
;The images here is twice the angular size so it's 0.05" per pixel
;cswa15 is at z=2.16 whose scale is 8.441 kpc/" or .422 kpc per pixel

outputdir='/scr2/nichal/workspace/paperplots/cswa15'   
;fitstoeps,file,scale,xtitle,ytitle,colorbartitle,directory_out
;Halpha
fitstoeps,'sourceha_pretty.fits',0.422,'kpc','cswa15','Ha flux(erg/s/cm2)',outputdir,format='(F6.3)'
;Halpha_detection
fitstoeps,'sourceha_detection_pretty.fits',0.422,'kpc','cswa15','Ha detection(sigma)',outputdir,format='(F6.3)'
;Nii
fitstoeps,'sourcenii_pretty.fits',0.422,'kpc','cswa15','[NII] flux(erg/s/cm2)',outputdir,format='(F6.3)'
;OIII
fitstoeps,'sourceoiii_pretty.fits',0.422,'kpc','cswa15','[OIII] flux(erg/s/cm2)',outputdir,format='(F6.3)'
;Hbeta
fitstoeps,'sourcehb_pretty.fits',0.422,'kpc','cswa15','Hb flux(erg/s/cm2)',outputdir,format='(F6.3)'
;kinematic
fitstoeps,'sourcekinematic_pretty.fits',0.422,'kpc','cswa15','velocity(km/s)',outputdir,format='(I0)'
;dispersion
fitstoeps,'sourceveldisp_pretty.fits',0.422,'kpc','cswa15','velocity dispersion(km/s)',outputdir,format='(I0)'
;N2index
fitstoeps,'sourcen2index_pretty.fits',0.422,'kpc','cswa15','[NII]/Ha',outputdir,format='(F5.2)'
;Bayesianmetal
;N2index
fitstoeps,'sourceBayesianmetal_pretty.fits',0.422,'kpc','cswa15','12+logO/H(M08)',outputdir,format='(F5.2)'

end
