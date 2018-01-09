pro add_psf,fwhm,cubefile,xc,yc
;add the measured psf to the Ha_cube
cube = readfits(cubefile,header)
sizec = size(cube)
ndim = min([sizec(1),sizec(2)])
model = psf_gaussian(npixel=ndim,ndimen=2,fwhm=fwhm,centroid=[xc,yc],/normalize)
psf = fltarr(sizec(1),sizec(2))
psf[0:ndim(0)-1,0:ndim(0)-1] = model

;overwrite if it exists. or add the the frame if it doesnt.
if sizec(3) ge 13 then cube[*,*,12] = psf else cube = [[[cube]],[[psf]]]
writefits,cubefile,cube,header

end
