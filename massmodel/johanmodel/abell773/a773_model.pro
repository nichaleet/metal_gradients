pro a773_model, image,hdr,pixsize,newimage,newheader,newscale
; IDL function to reconstruct the z=2.3 arc in Abell773 in the image plane.
;
; INPUTS:
; image   - a data cube with dimensions [x, y,*] ; note that usually osiris images have inaccurate RA-DEC, please use the already adjusted data or just use the extracted Ha,velocity map etc
; hdr     - header for a 2D slice (coordinates [x,y] only) of "image", with correct astrometric keywords
; pixsize - pixel size for "image" (arcseconds)
; OUTPUT:
; returns the source plane reconstruction, with a pixel scale of 0.005 arcseconds.


; xsky [ysky] is list of x [y] values of points on sky - these transform
; to position xsource [ysource] when the lensing potential is taken
; into account.
readcol, 'A773_source.dat',ids,xsource,ysource ;(id, RA,DEC)
readcol, 'A773_Osiris.mul',idi,xsky,ysky       ;image plane (idi, RA, DEC)

;The x,y given from Johan's model are too fine. Reduce to half the resolution. To be about 100 pixelx100 pixel which kinda match the osiris's resolution
ind_tobedelete = findgen(fix(max(ids)/1.5))*1.5
ind_tobedelete = fix(ind_tobedelete)
remove,ind_tobedelete,ids,xsource,ysource,idi,xsky,ysky


;WARNING: this lens model predicts multiple images for the source! Use only the data corresponding to the observed image,
; and ignore all other grid values.

;good = where(ysky lt 51.7263)
;idi = idi[good]
;xsky = xsky[good]
;ysky = ysky[good]
; Conveniently, "idi" and "ids" are now exactly the same!

; Convert xsky, ysky positions from (ra, dec) degrees into pixels
adxy, hdr, xsky, ysky, x, y
;stop
;now x,y is the pixel number of the coordinate given in image plane from johan's model

; Interpolate image from regular grid to xskypix,yskypix

; Interpolate image to points [x,y]
print, 'interpolating image'
;intimage = min_curve_surf(image, xgrid=[0,1], ygrid=[0,1], xpout=x, ypout=y)
; Can't use min_curve_surf because arrays are too big. Instead use
; bilinear interpolation which is very fast
intimage = interpolate(image[*,*,0], x, y,missing=1./0.)

; Map intimage values into source plane
   ; First triangulate the source-plane data points
triangulate, xsource, ysource, tr, b

   ; Determine the grid of x and y values to use in the source plane.
sourcescale = fltarr(n_elements(ysource)-1)
for i=0,n_elements(sourcescale)-1 do sourcescale(i)=ysource(i+1)-ysource(i)

;Choose one of the scales below by looking at the reasonable output
;pscale = ysource[1]-ysource[0]   ; pixel scale (arcseconds)
pscale = median(sourcescale)
pscale = 0.019

xgs = findgen(round((max(xsource)-min(xsource))/pscale) + 1.)*pscale + min(xsource)
ygs = findgen(round((max(ysource)-min(ysource))/pscale) + 1.)*pscale + min(ysource)
;xgs, ygs is the size of output x, y

   ; use "trigrid" to recover the 2D image
source_out = trigrid(xsource, ysource, intimage, tr, xout=xgs, yout=ygs)
help, xsource,ysource
help, xgs,ygs
;stop
; To conserve surface brightness,
; multiply by (source plane pixel size)/(image plane pixel size)
source_out = source_out * (pscale)^2. / (pixsize)^2.

; Now loop over all frames
source_cube = make_array((size(source_out))[1], (size(source_out))[2],(size(image))[3])
for i=0,(size(image))[3]-1 do begin
   intimage = interpolate(image[*,*,i], x, y)
   source_cube[*,*,i] = trigrid(xsource, ysource, intimage, tr, xout=xgs, yout=ygs)
endfor

print, 'scale is', pscale

mkhdr,newheader,source_cube
paraname = ['CRVAL1','CRVAL2','CRPIX1','CRPIX2','CDELT1','CDELT2']
paraval = [xsource[0],ysource[0],0,0,pscale/3600.,pscale/3600.]

hdr_paraname_str = ['CTYPE1','CTYPE2','CUNIT1','CUNIT2','RADESYS']
hdr_paraname_str_val = ['RA---TAN','DEC--TAN','deg','deg','FK5']

for ii=0,n_elements(hdr_paraname_Str)-1 do sxaddpar,newheader,hdr_paraname_str(ii),hdr_paraname_str_val(ii)
for ii=0,n_elements(paraname)-1 do sxaddpar, newheader,paraname[ii],paraval[ii]

newimage = source_cube
newscale= pscale
end
