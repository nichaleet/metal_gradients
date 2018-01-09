pro makeeps
;cswa28 is at z=2.09 whose scale is 8.478 kpc/" and the image used in model making is 0.1" per pixel. But the image is twice the angular size so the final scale is 0.424 kpc/pixel

outputdir='/scr2/nichal/workspace/paperplots/cswa28'   
;fitstoeps,file,scale,xtitle,ytitle,colorbartitle,directory_out

;Halpha
fitstoeps,'sourceha_pretty.fits',0.424,'kpc','cswa28','Ha flux(erg/s/cm2)',outputdir,format='(F6.3)'
;Halpha_detection
fitstoeps,'sourceha_detection_pretty.fits',0.424,'kpc','cswa28','Ha detection(sigma)',outputdir,format='(F6.3)'
;Nii
fitstoeps,'sourcenii_pretty.fits',0.424,'kpc','cswa28','[NII] flux(erg/s/cm2)',outputdir,format='(F6.3)'
;kinematic
fitstoeps,'sourcekinematic_pretty_corrected.fits',0.424,'kpc','cswa28','velocity(km/s)',outputdir,format='(I0)'
;dispersion
fitstoeps,'sourceveldisp_pretty.fits',0.424,'kpc','cswa28','velocity dispersion(km/s)',outputdir,format='(I0)'
;N2index
fitstoeps,'sourcen2index_pretty.fits',0.424,'kpc','cswa28','[NII]/Ha',outputdir,format='(F5.2)'

end
