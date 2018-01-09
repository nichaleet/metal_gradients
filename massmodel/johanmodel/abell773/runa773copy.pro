pro runa773copy
;
;;Re lense the hubble image
;imagehb = readfits('a773_small.fits',header)
;sizehb= size(imagehb)
;scale = 0.0995864
;imagehb= reform(imagehb,sizehb(1),sizehb(2),1)
;a773_model,imagehb,header,scale,newimage,newheader,newscale
;writefits,'a773_small_sourceplane.fits',newimage,newheader
;
;stop
;
;print, 'Now doing Halpha image '
;image =  readfits("/scr2/nichal/workspace/output/abell773_Ha_Kc3_mosaic_sky_230hr_acube.fits",header)
;scale = 0.1
;a773_model,image,header,scale,newimage,newheader,newscale
;
;newimage(where(newimage eq 0.)) = 1./0.
;writefits,"/scr2/nichal/workspace/output/abell773_Ha_Kc3_mosaic_sky_230hr_acube_sourceplane.fits",newimage,newheader
;
;;Make surface brightness conserved. Surface brightness is flux per solid angle
;; Now newimage(*,*,1) is flux in unit of erg/s/cm^2 so we have to devide the sum of flux by total angular area
;imbrightnew = newimage[*,*,1]*(newscale)^2./(scale)^2.
;
;;check surface brightness. Surface brightness is flux per solid angle
;; Now image(*,*,1) is flux in unit of erg/s/cm^2 so we have to devide the sum of flux by total angular area
;imbright = image(*,*,1)
;
;goodim = where(finite(imbright) eq 1.)
;goodimnew = where(finite(imbrightnew) eq 1.)
;brightness= total(imbright(goodim))/(n_elements(goodim)*scale^2)
;brightnessnew =  total(imbrightnew(goodimnew))/(n_elements(goodimnew)*newscale^2)
;print, brightness,brightnessnew
;
;
;print, 'Now doing Hbeta image '
;image =  readfits("/scr2/nichal/workspace/output/abell773_Hb_Hn3_mosaic_sky_2hr_acube_aligned.fits",header)
;scale = 0.1
;a773_model,image,header,scale,newimage,newheader,newscale
;
;
;newimage(where(newimage eq 0.)) = 1./0.
;writefits,"/scr2/nichal/workspace/output/abell773_Hb_Hn3_mosaic_sky_2hr_acube_aligned_sourceplane.fits",newimage,newheader
;
;;Make surface brightness conserved. Surface brightness is flux per solid angle
;; Now newimage(*,*,1) is flux in unit of erg/s/cm^2 so we have to devide the sum of flux by total angular area
;imbrightnew = newimage[*,*,1]*(newscale)^2./(scale)^2.
;
;;check surface brightness. Surface brightness is flux per solid angle
;; Now image(*,*,1) is flux in unit of erg/s/cm^2 so we have to devide the sum of flux by total angular area
;imbright = image(*,*,1)
;
;goodim = where(finite(imbright) eq 1.)
;goodimnew = where(finite(imbrightnew) eq 1.)
;brightness= total(imbright(goodim))/(n_elements(goodim)*scale^2)
;brightnessnew =  total(imbrightnew(goodimnew))/(n_elements(goodimnew)*newscale^2)
;print, brightness,brightnessnew
;
print, 'Now doing Metallicity image '
image =  readfits("/scr2/nichal/workspace/output/metallicity/Abell773_metallicity.fits",header)
scale = 0.1
a773_model,image,header,scale,newimage,newheader,newscale


newimage(where(newimage eq 0.)) = 1./0.
writefits,"/scr2/nichal/workspace/output/metallicity/Abell773_metallicity_sourceplane_finescale.fits",newimage,newheader

stop
end
;
