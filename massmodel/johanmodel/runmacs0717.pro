pro runmacs0717
;Re lense the hubble image
imagehb = readfits('macs0717_small.fits',header)
sizehb= size(imagehb)
scale = 0.065
imagehb= reform(imagehb,sizehb(1),sizehb(2),1) ;make it into a cube
macs0717_model,imagehb,header,scale,newimage,newheader,newscale
writefits,'macs0717_small_sourceplane.fits',newimage,newheader

print, 'Now doing Halpha image'
image =  readfits("/scr2/nichal/workspace/output/macs0717_Ha_Kc5_handmosaic_scaledsky_130hr_acube.fits",header)
;image(where(finite(image) eq 0.)) = 0.
scale = 0.1
macs0717_model,image,header,scale,newimage,newheader,newscale

newimage(where(newimage eq 0.)) = 1./0.
writefits,"/scr2/nichal/workspace/output/macs0717_Ha_Kc5_handmosaic_scaledsky_130hr_acube_sourceplane.fits",newimage,newheader

;Make surface brightness conserved. Surface brightness is flux per solid angle
; Now newimage(*,*,1) is flux in unit of erg/s/cm^2 so we have to devide the sum of flux by total angular area
imbrightnew = newimage[*,*,1]*(newscale)^2./(scale)^2.

;check surface brightness. Surface brightness is flux per solid angle
; Now image(*,*,1) is flux in unit of erg/s/cm^2 so we have to devide the sum of flux by total angular area
imbright = image(*,*,1)

goodim = where(finite(imbright) eq 1.)
goodimnew = where(finite(imbrightnew) eq 1.)
brightness= total(imbright(goodim))/(n_elements(goodim)*scale^2)
brightnessnew =  total(imbrightnew(goodimnew))/(n_elements(goodimnew)*newscale^2)
print, 'surface brightness before and after:'
print, brightness,brightnessnew

;print, 'Now doing Hbeta image'
;image =  readfits("/scr2/nichal/workspace/output/macs1133_Hb_Jn2_mosaic_1hr_acube.fits",header)
;scale = 0.1
;macs1133_model,image,header,scale,newimage,newheader,newscale

;newimage(where(newimage eq 0.)) = 1./0.
;writefits,"/scr2/nichal/workspace/output/macs1133_Hb_Jn2_mosaic_1hr_acube_sourceplane.fits",newimage,newheader

print, 'Now doing Metallicity image'
image =  readfits("/scr2/nichal/workspace/output/metallicity/macs0717_metallicitynew.fits",header)
scale = 0.1
macs0717_model,image,header,scale,newimage,newheader,newscale

newimage(where(newimage eq 0.)) = 1./0.
writefits,"/scr2/nichal/workspace/output/metallicity/macs0717_metallicitynew_sourceplane.fits",newimage,newheader


stop
end
