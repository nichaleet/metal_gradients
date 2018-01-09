pro cswa31_make_gal_region
;Read regions in Gemini
imref   = readfits('check_c31g.fits',hdr)
goodpix = where(imref eq 32,ctgood)
imsize = size(imref,/dimensions)
xarr = rebin(findgen(imsize(0)),imsize(0),imsize(1))
yarr = rebin(transpose(findgen(imsize(1))),imsize(0),imsize(1))

xgood = xarr(goodpix)
ygood = yarr(goodpix)

extast,hdr,astr
xy2ad,xgood,ygood,astr,ra,dec 

;convert to SDSS
im = readfits('cswa31sdss_u.fits',hdr)
extast,hdr,astr
imsize = size(im,/dimensions)
mask = fltarr(imsize(0),imsize(1))
ad2xy,ra,dec,astr,xout,yout
xout=round(xout)
yout=round(yout)

for i=0,ctgood-1 do mask[xout(i),yout(i)]=1.

writefits,'cswa31sdss_u_mask.fits',mask,hdr

stop
end
