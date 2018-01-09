pro oiii_psf_fix

im = readfits('sourceoiii_psf_pretty.fits')
im(where(finite(im) eq 0))=0.
;27,10-12
for ii=10,12 do im[27,ii] = total(im[26:28,ii-1:ii+1])/9.
for ii=10,12 do im[29,ii] = total(im[28:30,ii-1:ii+1])/9.
im(where(im eq 0.))= 1./0.
writefits,'sourceoiii_psf_pretty_fixed.fits',im
stop
end
