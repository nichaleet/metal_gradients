function undistort_cube_0744_singlet, image, hdr, pixsize

; IDL function to reconstruct the z=2.21 arc in MACS J0744 in the image plane.
;
; INPUTS:
; image   - a data cube with dimensions [lambda, x, y]
; hdr     - header for a 2D slice (coordinates [x,y] only) of "image", with correct astrometric keywords
; pixsize - pixel size for "image" (arcseconds)
; OUTPUT:
; returns the source plane reconstruction, with a pixel scale of 0.005 arcseconds.


; xsky [ysky] is list of x [y] values of points on sky - these transform
; to position xsource [ysource] when the lensing potential is taken
; into account.
readcol, '../lensmodels/macs0744_sourceinput.dat',ids,xsource,ysource
readcol,'../lensmodels/macs0744_sourceoutput.dat',idi,xsky, ysky

; WARNING: this lens model predicts multiple images for the source! Use only the data corresponding to the observed image,
; and ignore all other grid values.
good = where(ysky lt 39.46)
idi = idi[good]
xsky = xsky[good]
ysky = ysky[good]
; Conveniently, "idi" and "ids" are now exactly the same!

; Convert xsky, ysky positions from (ra, dec) degrees into pixels
adxy, hdr, xsky, ysky, x, y

; Interpolate image from regular grid to xskypix,yskypix

; Interpolate image to points [x,y]
print, 'interpolating image'
;intimage = min_curve_surf(image, xgrid=[0,1], ygrid=[0,1], xpout=x, ypout=y)
; Can't use min_curve_surf because arrays are too big. Instead use
; bilinear interpolation which is very fast
intimage = interpolate(image[0,*,*], x, y)

; Map intimage values into source plane
   ; First triangulate the source-plane data points
triangulate, xsource, ysource, tr, b
   ; Determine the grid of x and y values to use in the source plane.
pscale = ysource[1]-ysource[0]   ; pixel scale (arcseconds)
xgs = findgen(round((max(xsource)-min(xsource))/pscale) + 1.)*pscale + min(xsource)
ygs = findgen(round((max(ysource)-min(ysource))/pscale) + 1.)*pscale + min(ysource)
   ; use "trigrid" to recover the 2D image
source_out = trigrid(xsource, ysource, intimage, tr, xout=xgs, yout=ygs)

; To conserve surface brightness,
; multiply by (source plane pixel size)/(image plane pixel size)
source_out = source_out * (pscale)^2. / (pixsize)^2.

; Now loop over all wavelengths
source_cube = make_array((size(image))[1], (size(source_out))[1], (size(source_out))[2])
for i=0,(size(image))[1]-1 do begin
   intimage = interpolate(image[i,*,*], x, y)
   source_cube[i,*,*] = trigrid(xsource, ysource, intimage, tr, xout=xgs, yout=ygs)
endfor


return, source_cube

end
