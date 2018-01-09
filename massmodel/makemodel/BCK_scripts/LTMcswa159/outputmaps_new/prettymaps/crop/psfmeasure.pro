pro psfmeasure
im=readfits('sourceha_psf_pretty.fits')
im(where(~finite(im)))=0.
im=im[8:28,14:24]
;smooth
im = smooth(im,2)
results=gauss2dfit(im,param,/tilt)

widthx=sqrt(param(2)^2)
widthy=sqrt(param(3)^2)
print,'psf width', [param(2),param(3)]*0.418205

writefits,'psf_chopped.fits',[[[im]],[[results]]]
stop
end
