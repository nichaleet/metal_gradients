pro makeeps
;cswa31  the final scale is 0.628 kpc/pixel

outputdir='/scr2/nichal/workspace/paperplots/cswa31'   
;fitstoeps,file,scale,xtitle,ytitle,colorbartitle,directory_out

;Halpha
fitstoeps,'sourceha_pretty.fits',0.628,'kpc','cswa31','Ha flux(erg/s/cm2)',outputdir,format='(F6.3)'
;Halpha_detection
fitstoeps,'sourceha_detection_pretty.fits',0.628,'kpc','cswa31','Ha detection(sigma)',outputdir,format='(F6.3)'
;Nii
fitstoeps,'sourcenii_pretty.fits',0.628,'kpc','cswa31','[NII] flux(erg/s/cm2)',outputdir,format='(F6.3)'
;kinematic
fitstoeps,'sourcekinematic_pretty_corrected.fits',0.628,'kpc','cswa31','velocity(km/s)',outputdir,format='(I0)'
;dispersion
fitstoeps,'sourceveldisp_pretty.fits',0.628,'kpc','cswa31','velocity dispersion(km/s)',outputdir,format='(I0)'
;N2index
fitstoeps,'sourcen2index_pretty.fits',0.628,'kpc','cswa31','[NII]/Ha',outputdir,format='(F5.2)'

end
