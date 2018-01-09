pro abell773_make_gal_region
;Read regions in HST
imref   = mrdfits('check_abell_hst.fits',1,hdr)
goodpix = where(imref eq 524,ctgood)
imsize = size(imref,/dimensions)
xarr = rebin(findgen(imsize(0)),imsize(0),imsize(1))
yarr = rebin(transpose(findgen(imsize(1))),imsize(0),imsize(1))

xgood = xarr(goodpix)
ygood = yarr(goodpix)

extast,hdr,astr
xy2ad,xgood,ygood,astr,ra,dec 

;convert to SDSS
im = readfits('a773sdss_g.fits',hdr)
extast,hdr,astr
imsize = size(im,/dimensions)
mask = fltarr(imsize(0),imsize(1))
ad2xy,ra,dec,astr,xout,yout
xout=round(xout)
yout=round(yout)

for i=0,ctgood-1 do mask[xout(i),yout(i)]=1.

writefits,'a773sdss_g_mask.fits',mask,hdr

stop
end
